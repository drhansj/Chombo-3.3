#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "InterfaceJump.H"
#include "EBArith.H"
#include "BaseIVFactory.H"
#include "EBLevelGrid.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"

InterfaceJump::InterfaceJump()
{
  m_isDefined = false;
}

InterfaceJump::~InterfaceJump()
{
}

void InterfaceJump::define(const MFIndexSpace      & a_MFIS,
                           const ProblemDomain     & a_domain,
                           const DisjointBoxLayout & a_dbl,
                           const RealVect          & a_vectDx,
                           const RealVect          & a_physOrigin,
                           const bool                a_isVectorJump)
{
  CH_TIME("InterfaceJump::define");
  m_dbl = a_dbl;
  m_domain = a_domain;
  m_vectDx = a_vectDx;
  int numEBISLGhost = 3;
  m_nComp = 1;
  m_origin = a_physOrigin;
  m_isVectorJump = a_isVectorJump;

  //mike's anisotropic dx
  setDxConstants();

  int phase = 0;
  definePhase(a_MFIS,
              a_domain,
              m_dbl,
              m_0ebisl,
              m_0stencils,
              m_0inHomDirichWeight,
              m_0dirichletDropOrder,
              m_0boundaryVofs,
              m_0boundaryFaceIndices,
              numEBISLGhost,
              phase);

  phase = 1;
  definePhase(a_MFIS,
              a_domain,
              m_dbl,
              m_1ebisl,
              m_1stencils,
              m_1inHomDirichWeight,
              m_1dirichletDropOrder,
              m_1boundaryVofs,
              m_1boundaryFaceIndices,
              numEBISLGhost,
              phase);

  if (a_isVectorJump)
    {
      defineVectorJumpLayoutData(m_0ebisl,
                                 m_vectorGD,
                                 m_vectorGN);
    }
  else
    {
      defineScalarJumpLayoutData(m_0ebisl,
                                 m_scalarGD,
                                 m_scalarGN);
    }

  computeBoundaryLayoutData(m_0boundaryVofs,
                            m_1boundaryVofs,
                            m_0boundaryFaceIndices,
                            m_1boundaryFaceIndices);

  cacheVoFIt();
}

void InterfaceJump::
defineScalarJumpLayoutData(EBISLayout&                               a_ebisl,
                           LayoutData< BaseIVFAB< Vector<Real> > >&  a_scalarGD,
                           LayoutData< BaseIVFAB< Vector<Real> > >&  a_scalarGN)
{
  CH_TIME("InterfaceJump::defineScalarJumpLayoutData");
  a_scalarGD.define(m_dbl);
  a_scalarGN.define(m_dbl);
  for (DataIterator dit = m_dbl.dataIterator(); dit.ok(); ++dit)
    {
      Box              region      = m_dbl.get(dit());
      const EBISBox    ebisBox     = a_ebisl[dit()];
      IntVectSet       ivsIrreg    = ebisBox.getIrregIVS(region);
      IntVectSet       ivsMulti    = ebisBox.getMultiCells(region);

      a_scalarGD[dit()].define(ivsIrreg, ebisBox.getEBGraph(), 1);
      a_scalarGN[dit()].define(ivsIrreg, ebisBox.getEBGraph(), 1);

      for (VoFIterator vofIt(ivsIrreg,ebisBox.getEBGraph()); vofIt.ok(); ++vofIt)
        {
          VolIndex vof = vofIt();
          int nFace = ebisBox.numFacePhase(vof);
          a_scalarGD[dit()](vof, 0).resize(nFace, 0.0);
          a_scalarGN[dit()](vof, 0).resize(nFace, 0.0);
        }
    }
}

void InterfaceJump::
defineVectorJumpLayoutData(EBISLayout&                                   a_ebisl,
                           LayoutData< BaseIVFAB< Vector<RealVect> > >&  a_vectorGD,
                           LayoutData< BaseIVFAB< Vector<RealVect> > >&  a_vectorGN)
{
  CH_TIME("InterfaceJump::defineVectorJumpLayoutData");
  a_vectorGD.define(m_dbl);
  a_vectorGN.define(m_dbl);
  for (DataIterator dit = m_dbl.dataIterator(); dit.ok(); ++dit)
    {
      Box              region      = m_dbl.get(dit());
      const EBISBox    ebisBox     = a_ebisl[dit()];
      IntVectSet       ivsIrreg    = ebisBox.getIrregIVS(region);
      IntVectSet       ivsMulti    = ebisBox.getMultiCells(region);

      a_vectorGD[dit()].define(ivsIrreg, ebisBox.getEBGraph(), 1);
      a_vectorGN[dit()].define(ivsIrreg, ebisBox.getEBGraph(), 1);

      for (VoFIterator vofIt(ivsIrreg,ebisBox.getEBGraph()); vofIt.ok(); ++vofIt)
        {
          VolIndex vof = vofIt();
          int nFace = ebisBox.numFacePhase(vof);
          a_vectorGD[dit()](vof, 0).resize(nFace, RealVect::Zero);
          a_vectorGN[dit()](vof, 0).resize(nFace, RealVect::Zero);
        }
    }
}

void InterfaceJump::definePhase(const MFIndexSpace&                             a_MFIS,
                                const ProblemDomain&                            a_domain,
                                const DisjointBoxLayout&                        a_dbl,
                                EBISLayout&                                     a_ebisl,
                                LayoutData< BaseIVFAB< Vector<VoFStencil> > >&  a_stencils,
                                LayoutData< BaseIVFAB< Vector<Real> > >&        a_inHomDirichWeight,
                                LayoutData< BaseIVFAB< Vector<int> > >&         a_dirichletDropOrder,
                                LayoutData< BaseIVFAB< Vector<VolIndex> > >&    a_boundaryVofs,
                                LayoutData< BaseIVFAB< Vector<int> > >&         a_boundaryFaceIndices,
                                const int&                                      a_numEBISLGhost,
                                const int&                                      a_phase)
{
  CH_TIME("InterfaceJump::definePhase");
  a_MFIS.fillEBISLayout(a_ebisl,
                        a_phase,
                        m_dbl,
                        a_domain.domainBox(),
                        a_numEBISLGhost);

  defineLayoutData(a_ebisl,
                   a_stencils,
                   a_inHomDirichWeight,
                   a_dirichletDropOrder,
                   a_boundaryVofs,
                   a_boundaryFaceIndices);

  for (DataIterator dit = m_dbl.dataIterator(); dit.ok(); ++dit)
    {
      Box  region  = m_dbl.get(dit());
      const EBISBox ebisBox = a_ebisl[dit()];

      //one irregular ivs for both phases
      IntVectSet ivsIrreg  = ebisBox.getIrregIVS(region);

      BaseIVFAB< Vector<Real> >&        inHomWeightFab = a_inHomDirichWeight[dit()];
      BaseIVFAB< Vector<int> >&         dropOrder = a_dirichletDropOrder[dit()];
      BaseIVFAB< Vector<VoFStencil> >&  ivStencil = a_stencils[dit()];

      for (VoFIterator ivof(ivsIrreg,ebisBox.getEBGraph()); ivof.ok(); ++ivof)
        {
          VolIndex vof = ivof();

          int nFace = ebisBox.numFacePhase(vof);
          inHomWeightFab(vof, 0).resize(nFace);
          dropOrder(vof, 0).resize(nFace);
          ivStencil(vof, 0).resize(nFace);

          for (int iface=0; iface<nFace; iface++)
            {
              VoFStencil& stencil = ivStencil(vof,0)[iface];
              Real bdArea = ebisBox.bndryArea(vof, iface);
              //if not dirichletDropOrder then stenSize = 2*SpaceDim+1
              int dirichletDropOrder = 999;

              //this is for when the EB is Dirichlet
              Real inHomDirich       = 999.999999;

              VoFStencil johanStencil;
              RealVect normal = ebisBox.normal(vof, iface);
              RealVect centroid = ebisBox.bndryCentroid(vof, iface);
              johanDirichletStencil(normal,
                                    centroid,
                                    vof,
                                    ebisBox,
                                    johanStencil,
                                    dirichletDropOrder,
                                    inHomDirich);

              if (dirichletDropOrder == -1)
                {
                  stencil += johanStencil;
                }
              else if (dirichletDropOrder == 1)
                {
                  int ivar = 0;
                  VoFStencil lsStencil;
                  leastSquaresDirichStencil(lsStencil,
                                            inHomDirich,
                                            normal,
                                            centroid,
                                            bdArea,
                                            vof,
                                            ebisBox,
                                            m_vectDx,
                                            m_domain,
                                            ivar);
                  if (lsStencil.size() > 0)
                    {
                      dirichletDropOrder = -1;
                      stencil += lsStencil;
                    }
                  else
                    {
#ifdef NDEBUG
                      MayDay::Warning("InterfaceJump::definePhase - empty stencil");
#else
                      if (bdArea > 0)
                        {
                          char msg[1024];
                          IntVect iv = vof.gridIndex();
                          IntVect smallEnd(a_domain.domainBox().smallEnd());
                          IntVect   bigEnd(a_domain.domainBox().bigEnd());
                          if (SpaceDim == 2)
                            {
                              sprintf(msg,
                                      "InterfaceJump::definePhase - empty stencil in phase %d at (%d,%d) on domain box (%d,%d) to (%d,%d)",
                                      a_phase,iv[0],iv[1],
                                      smallEnd[0],smallEnd[1],bigEnd[0],bigEnd[1]);
                            }
                          else if (SpaceDim == 3)
                            {
                              sprintf(msg,
                                      "InterfaceJump::definePhase - empty stencil in phase %d at (%d,%d,%d) on domain box (%d,%d,%d) to (%d,%d,%d)",
                                      a_phase,iv[0],iv[1],iv[2],
                                      smallEnd[0],smallEnd[1],smallEnd[2],
                                      bigEnd[0],bigEnd[1],bigEnd[2]);
                            }
                          else
                            {
                              sprintf(msg,"InterfaceJump::definePhase - empty stencil in phase %d",
                                      a_phase);
                            }
                          MayDay::Warning(msg);
                        }
#endif
                    }
                }

              //this = -1 unless we can't even do least squares
              dropOrder(vof,0)[iface] = dirichletDropOrder;
              //this will be call zeta later
              inHomWeightFab(vof,0)[iface] = inHomDirich;
              if (bdArea > 0.0)
                {
                  stencil *= m_vectDx[0]/bdArea;
                  inHomWeightFab(vof,0)[iface] *= m_vectDx[0]/bdArea;
                }
              else
                {
#ifdef NDEBUG
                  MayDay::Warning("InterfaceJump::definePhase - boundary with zero boundary area");
#else
                  char msg[1024];
                  IntVect iv = vof.gridIndex();
                  IntVect smallEnd(a_domain.domainBox().smallEnd());
                  IntVect   bigEnd(a_domain.domainBox().bigEnd());
                  if (SpaceDim == 2)
                    {
                      sprintf(msg,
                              "InterfaceJump::definePhase - zero boundary area in phase %d at (%d,%d) on domain box (%d,%d) to (%d,%d)",
                              a_phase,iv[0],iv[1],
                              smallEnd[0],smallEnd[1],bigEnd[0],bigEnd[1]);
                    }
                  else if (SpaceDim == 3)
                    {
                      sprintf(msg,
                              "InterfaceJump::definePhase - zero boundary area in phase %d at (%d,%d,%d) on domain box (%d,%d,%d) to (%d,%d,%d)",
                              a_phase,iv[0],iv[1],iv[2],
                              smallEnd[0],smallEnd[1],smallEnd[2],
                              bigEnd[0],bigEnd[1],bigEnd[2]);
                    }
                  else
                    {
                      sprintf(msg,"InterfaceJump::definePhase - zero boundary area in phase %d",
                              a_phase);
                    }
                  MayDay::Warning(msg);
#endif
                }
            }
        }
    }
  m_isDefined = true;
}

void InterfaceJump::leastSquaresDirichStencil(VoFStencil&          a_stencil,
                                              Real&                a_weight,
                                              const RealVect&      a_normal,
                                              const RealVect&      a_centroid,
                                              const Real&          a_bndryArea,
                                              const VolIndex&      a_vof,
                                              const EBISBox&       a_ebisBox,
                                              const RealVect&      a_dx,
                                              const ProblemDomain& a_domain,
                                              int                  a_ivar)
{
  CH_TIME("InterfaceJump::leastSquaresDirichStencil");
  Real scaling = m_dxScale*a_bndryArea;

  EBArith::getLeastSquaresGradSten(a_stencil, a_weight, a_normal, a_centroid,
                                   a_vof, a_ebisBox, a_dx, a_domain, a_ivar);

  if (a_stencil.size() == 0)
    {
#ifdef NDEBUG
      MayDay::Warning("InterfaceJump::definePhase - trying all VoFs stencil.");
#else
      char msg[1024];
      IntVect iv = a_vof.gridIndex();
      IntVect smallEnd(a_domain.domainBox().smallEnd());
      IntVect   bigEnd(a_domain.domainBox().bigEnd());
      if (SpaceDim == 2)
        {
          sprintf(msg,
                  "InterfaceJump::definePhase - trying all VoFs stencil at (%d,%d) on domain box (%d,%d) to (%d,%d)",
                  iv[0],iv[1],
                  smallEnd[0],smallEnd[1],bigEnd[0],bigEnd[1]);
        }
      else if (SpaceDim == 3)
        {
          sprintf(msg,
                  "InterfaceJump::definePhase - trying all VoFs stencil at (%d,%d,%d) on domain box (%d,%d,%d) to (%d,%d,%d)",
                  iv[0],iv[1],iv[2],
                  smallEnd[0],smallEnd[1],smallEnd[2],
                  bigEnd[0],bigEnd[1],bigEnd[2]);
        }
      else
        {
          sprintf(msg,"InterfaceJump::definePhase - trying all VoFs stencil");
        }
      MayDay::Warning(msg);
#endif
      EBArith::getLeastSquaresGradStenAllVoFs(a_stencil, a_weight, a_normal,
                                              a_centroid, a_vof, a_ebisBox,
                                              a_dx, a_domain, a_ivar);
    }

  a_stencil *= scaling;
  a_weight  *= -scaling;
}

void InterfaceJump::johanDirichletStencil(const RealVect& a_normal,
                                          const RealVect& a_centroid,
                                          const VolIndex& a_vof,
                                          const EBISBox&  a_ebisBox,
                                          VoFStencil&     a_vofStencil,
                                          int&            a_dirichletDropOrder,
                                          Real&           a_weight)
{
  CH_TIME("InterfaceJump::johanDirichStencil");
  Vector<VoFStencil> pointStencils;
  Vector<Real>       distanceAlongLine;
  IntVectSet cfivs;
  bool dropOrder = false;
  EBArith::johanStencil(dropOrder, pointStencils, distanceAlongLine,
                        a_normal, a_centroid, a_vof,
                        a_ebisBox, m_vectDx, cfivs);

  a_vofStencil.clear();
  if (dropOrder)
    {
      a_dirichletDropOrder = 1;
      return;
    }
  else
    {
      //if we got this far, sizes should be at least 2
      CH_assert(distanceAlongLine.size() >= 2);
      CH_assert(pointStencils.size() >= 2);
      Real x1 = distanceAlongLine[0];
      Real x2 = distanceAlongLine[1];
      //fit quadratic function to line and find the gradient at the origin
      //grad = (x2*x2*(phi1-phi0)-x1*x1(phi2-phi0))/(x2*x2*x1 - x1*x1*x2);
      Real denom = x2*x2*x1 - x1*x1*x2;
      //not done by reference because i want point stencils checkable externally.
      VoFStencil phi1Sten = pointStencils[0];
      VoFStencil phi2Sten = pointStencils[1];
      Real scaling = m_dxScale*a_ebisBox.bndryArea(a_vof);
      phi1Sten *=-scaling*x2*x2/denom;
      phi2Sten *= scaling*x1*x1/denom;
      //weight is the multiplier of the inhomogeneous value (phi0)
      a_weight =-scaling*(-x1*x1/denom + x2*x2/denom);
      a_vofStencil += phi1Sten;
      a_vofStencil += phi2Sten;
      a_dirichletDropOrder = -1;
    }
}

void InterfaceJump::boundarySlope(Real*            a_dPhiDn,
                                  const RealVect&  a_normal,
                                  const VolIndex&  a_vof,
                                  const EBCellFAB& a_phiFab,
                                  const EBISBox&   a_ebisBox,
                                  const int&       a_ivar)
{
  CH_TIME("InterfaceJump::boundarySlope");
  a_dPhiDn[0]=0;
  RealVect grad;
  for (int i=0; i<CH_SPACEDIM; ++i)
    {
      IntVect iv = a_vof.gridIndex();

      iv[i]+=1;

      int dir=1;
      int k = (i+1) % CH_SPACEDIM;
      if (!m_domain.contains(iv) || a_ebisBox.isCovered(iv))
      {
        iv[i]-=2;
        dir=-1;
      }
      if (!m_domain.contains(iv) ||a_ebisBox.isCovered(iv))
      {
        iv[k]+=1;
      }
      if (!m_domain.contains(iv) ||a_ebisBox.isCovered(iv))
      {
        iv[k]-=2;
      }
      if (!m_domain.contains(iv) ||a_ebisBox.isCovered(iv))
      {
        iv[i]+=2;
        dir=1;
      }
      if (!m_domain.contains(iv) ||a_ebisBox.isCovered(iv))
      {
        iv[k]+=2;
      }
      if (!m_domain.contains(iv) ||a_ebisBox.isCovered(iv))
      {
        iv[k]-=1;
        k=(i+2) % CH_SPACEDIM;
        iv[k]+=1;
      }
      if (!m_domain.contains(iv) ||a_ebisBox.isCovered(iv))
      {
        iv[k]-=2;
      }
      if (!m_domain.contains(iv) ||a_ebisBox.isCovered(iv))
      {
        iv[i]-=2;
        dir=-1;
      }

      CH_assert(!a_ebisBox.isCovered(iv));
      grad[i] = dir*(a_phiFab(VolIndex(iv,0),a_ivar) - a_phiFab(a_vof, a_ivar))/m_vectDx[i];
      a_dPhiDn[0] += a_normal[i]*grad[i];
    }
}

void InterfaceJump::
computeBoundaryLayoutData(LayoutData< BaseIVFAB< Vector<VolIndex> > >&  a_0boundaryVofs,
                          LayoutData< BaseIVFAB< Vector<VolIndex> > >&  a_1boundaryVofs,
                          LayoutData< BaseIVFAB< Vector<int> > >&       a_0boundaryFaceIndices,
                          LayoutData< BaseIVFAB< Vector<int> > >&       a_1boundaryFaceIndices)
{
  CH_TIME("InterfaceJump::computeBoundaryLayoutData");
  int nPhase = 2;

  // need EBISBox for boundary areas
  Vector<EBISBox> ebisBoxV;
  ebisBoxV.resize(nPhase);

  Vector< BaseIVFAB< Vector<VolIndex> >* > boundaryVofFabs(nPhase);
  Vector< BaseIVFAB< Vector<int> >* > boundaryFaceIndexFabs(nPhase);

  for (DataIterator dit = m_dbl.dataIterator(); dit.ok(); ++dit)
    {
      ebisBoxV[0] = m_0ebisl[dit()];
      ebisBoxV[1] = m_1ebisl[dit()];

      Vector<VoFIterator> vofItV;
      vofItV.resize(nPhase);
      for (int iphase=0; iphase<nPhase; iphase++)
        {
          Box         region    = m_dbl.get(dit());
          IntVectSet  ivsBndry  = ebisBoxV[iphase].boundaryIVS(region);
          EBGraph     ebGraph   = ebisBoxV[iphase].getEBGraph();
          vofItV[iphase].define(ivsBndry, ebGraph);
        }

      boundaryVofFabs[0] = &(a_0boundaryVofs[dit()]);
      boundaryVofFabs[1] = &(a_1boundaryVofs[dit()]);
      boundaryFaceIndexFabs[0] = &(a_0boundaryFaceIndices[dit()]);
      boundaryFaceIndexFabs[1] = &(a_1boundaryFaceIndices[dit()]);

      for (int iphase=0; iphase<nPhase; iphase++)
        {
          // index of the other phase; assume two phases
          int jphase = 1-iphase;

          for (vofItV[iphase].reset(); vofItV[iphase].ok();
               ++vofItV[iphase])
            {
              VolIndex ivof = vofItV[iphase]();
              int niFace = ebisBoxV[iphase].numFacePhase(ivof);

              // get the VolIndexes for the VoFs on my face(s)
              for (int iface=0; iface<niFace; iface++)
                {
                  VolIndex jvof = ebisBoxV[iphase].faceIndex(ivof, iface);
                  int jface = -1;
                  for (int jjface=0; jjface<ebisBoxV[jphase].numFacePhase(jvof); jjface++)
                    {
                      if (ebisBoxV[jphase].faceIndex(jvof, jjface) == ivof)
                        {
                          jface = jjface;
                        }
                    }
                  if (jface == -1)
                    {
                      MayDay::Error(" InterfaceJump::computeBoundaryLayoutData - unable to find matching face");
                    }
                  (*boundaryVofFabs[iphase])(ivof, 0)[iface] = jvof;
                  (*boundaryFaceIndexFabs[iphase])(ivof, 0)[iface] = jface;
                }
            }
        }
    }
}

void InterfaceJump::
defineLayoutData(EBISLayout&                                     a_ebisl,
                 LayoutData< BaseIVFAB< Vector<VoFStencil> > >&  a_stencils,
                 LayoutData< BaseIVFAB< Vector<Real> > >&        a_inHomDirichWeight,
                 LayoutData< BaseIVFAB< Vector<int> > >&         a_dirichletDropOrder,
                 LayoutData< BaseIVFAB< Vector<VolIndex> > >&    a_boundaryVofs,
                 LayoutData< BaseIVFAB< Vector<int> > >&         a_boundaryFaceIndices)
{
  CH_TIME("InterfaceJump::defineLayoutData");
  a_dirichletDropOrder.define(m_dbl);
  a_stencils.define(m_dbl);
  a_inHomDirichWeight.define(m_dbl);
  a_boundaryVofs.define(m_dbl);
  a_boundaryFaceIndices.define(m_dbl);

  m_vofItBndry.define(m_dbl);

  for (DataIterator dit = m_dbl.dataIterator(); dit.ok(); ++dit)
    {
      Box              region      = m_dbl.get(dit());
      const EBISBox    ebisBox     = a_ebisl[dit()];
      IntVectSet       ivsIrreg    = ebisBox.getIrregIVS(region);
      IntVectSet       ivsMulti    = ebisBox.getMultiCells(region);

      // Gets the cells for which it is nedeed to drop order
      a_dirichletDropOrder[dit()].define(ivsIrreg,ebisBox.getEBGraph(),1);

      // Define/allocate the BaseIVFab< VoFStencil >
      a_stencils         [dit()].define(ivsIrreg,ebisBox.getEBGraph(),1);
      a_inHomDirichWeight[dit()].define(ivsIrreg,ebisBox.getEBGraph(),1);

      // Define the Vector<VolIndex>s & Vector<int>s for boundary calculations
      a_boundaryVofs       [dit()].define(ivsIrreg,ebisBox.getEBGraph(),1);
      a_boundaryFaceIndices[dit()].define(ivsIrreg,ebisBox.getEBGraph(),1);

      for (VoFIterator vofIt(ivsIrreg,ebisBox.getEBGraph()); vofIt.ok(); ++vofIt)
        {
          VolIndex vof = vofIt();
          int nFace = ebisBox.numFacePhase(vof);
          a_stencils[dit()](vof, 0).resize(nFace);
          a_inHomDirichWeight[dit()](vof, 0).resize(nFace);
          a_dirichletDropOrder[dit()](vof, 0).resize(nFace);
          a_boundaryVofs[dit()](vof, 0).resize(nFace);
          a_boundaryFaceIndices[dit()](vof, 0).resize(nFace);
        }
    }
}

void InterfaceJump::cacheVoFIt()
{
  CH_TIME("InterfaceJump::cacheVoFIt");

  m_vofItBndry.define(m_dbl);

  for (DataIterator dit = m_dbl.dataIterator(); dit.ok(); ++dit)
    {
      Box              region      = m_dbl.get(dit());
      const EBISBox    ebisBox     = m_0ebisl[dit()];
      IntVectSet       ivsBndry    = ebisBox.boundaryIVS(region);
      //cache the vofIterators
      m_vofItBndry[dit()].define(ivsBndry,ebisBox.getEBGraph());
    }
}

void InterfaceJump::setDxConstants()
{
  //pre-compute (dx*dx, dy*dy,...)
  for (int idir=0;idir<SpaceDim;idir++)
    {
      Real power = 2.0;
      m_d2Vect[idir] = pow(m_vectDx[idir],power);
    }

  //compute m_dxScale (used for unscaling the boundary area)
  Real maxDx = 0.0;
  Real dv    = 1.0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      dv *= m_vectDx[idir];
      if ( m_vectDx[idir] > maxDx )
        {
          maxDx = m_vectDx[idir];
        }
    }
  m_dxScale = pow(maxDx,SpaceDim-1)/dv;

}

void InterfaceJump::setJump(const Real& a_gD,
                            const Real& a_gN)
{
  CH_TIME("InterfaceJump::setJump(const Real)");
  CH_assert(!m_isVectorJump);
  for (DataIterator dit = m_dbl.dataIterator(); dit.ok(); ++dit)
    {
      BaseIVFAB< Vector<Real> >& dirFab = m_scalarGD[dit()];
      BaseIVFAB< Vector<Real> >& neuFab = m_scalarGN[dit()];

      const Box region = m_dbl.get(dit());
      const EBISBox ebisBox = m_0ebisl[dit()];
      IntVectSet ivsIrreg = ebisBox.getIrregIVS(region);

      for (VoFIterator vofIt(ivsIrreg, ebisBox.getEBGraph()); vofIt.ok(); ++vofIt)
        {
          VolIndex vof = vofIt();
          int nFace = ebisBox.numFacePhase(vof);
          for (int iface=0; iface<nFace; iface++)
            {
              dirFab(vof, 0)[iface] = a_gD;
              neuFab(vof, 0)[iface] = a_gN;
            }
        }
    }
}

void InterfaceJump::setJump(const RealVect& a_gD,
                            const RealVect& a_gN)
{
  CH_TIME("InterfaceJump::setJump(const RealVect)");
  CH_assert(m_isVectorJump);
  for (DataIterator dit = m_dbl.dataIterator(); dit.ok(); ++dit)
    {
      BaseIVFAB< Vector<RealVect> >& dirFab = m_vectorGD[dit()];
      BaseIVFAB< Vector<RealVect> >& neuFab = m_vectorGN[dit()];

      const Box region = m_dbl.get(dit());
      const EBISBox ebisBox = m_0ebisl[dit()];
      IntVectSet ivsIrreg = ebisBox.getIrregIVS(region);

      for (VoFIterator vofIt(ivsIrreg, ebisBox.getEBGraph()); vofIt.ok(); ++vofIt)
        {
          VolIndex vof = vofIt();
          int nFace = ebisBox.numFacePhase(vof);
          for (int iface=0; iface<nFace; iface++)
            {
              dirFab(vof, 0)[iface] = a_gD;
              neuFab(vof, 0)[iface] = a_gN;
            }
        }
    }
}

void InterfaceJump::
setJump(const Vector< RefCountedPtr<BaseBCValue> >& a_phiValVect,
        const Vector< RefCountedPtr<BaseBCValue> >& a_flxValVect,
        const Vector<Real>&                         a_alpha,
        const Vector<Real>&                         a_beta,
        const Real&                                 a_time)
{
  CH_TIME("InterfaceJump::setJump(BaseBCValue)");
  CH_assert(!m_isVectorJump);
  CH_assert((a_alpha.size()==2) && (a_beta.size()==2));

  m_phiValVect = a_phiValVect;
  m_flxValVect = a_flxValVect;

  calculateJump(a_alpha, a_beta, a_time);
}

void InterfaceJump::
calculateJump(const Vector<Real>&  a_alpha,
              const Vector<Real>&  a_beta,
              const Real&          a_time)
{
  CH_TIME("InterfaceJump::calculateJump");
  for (DataIterator dit = m_dbl.dataIterator(); dit.ok(); ++dit)
    {
      BaseIVFAB< Vector<Real> >& dirFab = m_scalarGD[dit()];
      BaseIVFAB< Vector<Real> >& neuFab = m_scalarGN[dit()];

      const Box region = m_dbl.get(dit());
      const EBISBox ebisBox = m_0ebisl[dit()];
      IntVectSet ivsIrreg = ebisBox.getIrregIVS(region);

      for (VoFIterator vofIt(ivsIrreg, ebisBox.getEBGraph()); vofIt.ok(); ++vofIt)
        {
          VolIndex vof = vofIt();
          int nFace = ebisBox.numFacePhase(vof);
          for (int iface=0; iface<nFace; iface++)
            {
              RealVect normal = ebisBox.normal(vof, iface);
              RealVect loc = vof.gridIndex();
              loc += 0.5;
              loc += ebisBox.bndryCentroid(vof, iface);
              loc *= m_vectDx;
              loc += m_origin;
              dirFab(vof, 0)[iface] = m_phiValVect[1]->value(loc, RealVect::Zero, a_time, 0)
                                    - m_phiValVect[0]->value(loc, RealVect::Zero, a_time, 0);
              neuFab(vof, 0)[iface] = a_beta[1]*m_flxValVect[1]->value(loc, normal, a_time, 0)
                                    - a_beta[0]*m_flxValVect[0]->value(loc, normal, a_time, 0);
            }
        }
    }
}

void InterfaceJump::resetJump(const Vector<Real>&  a_alpha,
                              const Vector<Real>&  a_beta,
                              const Real&          a_time)
{
  if ( (m_phiValVect.size() == 2) &&
       (m_flxValVect.size() == 2) )
    {
      calculateJump(a_alpha, a_beta, a_time);
    }
}

void InterfaceJump::getFaceStencil(VoFStencil&          a_istencil,
                                   VoFStencil&          a_jstencil,
                                   const Vector<Real>&  a_beta,
                                   const VolIndex&      a_ivof,
                                   const DataIndex&     a_dix,
                                   const int&           a_iphase,
                                   const int&           a_iface)
{
  VolIndex jvof;
  int jface;
  Real iweight, jweight;
  int idrop, jdrop;
  a_istencil.clear();
  a_jstencil.clear();
  if (a_iphase == 0)
    {
      jvof = m_0boundaryVofs[a_dix](a_ivof,0)[a_iface];
      jface = m_0boundaryFaceIndices[a_dix](a_ivof,0)[a_iface];

      a_istencil += m_0stencils[a_dix](a_ivof, 0)[a_iface];
      a_jstencil += m_1stencils[a_dix](  jvof, 0)[  jface];

      iweight = m_0inHomDirichWeight[a_dix](a_ivof, 0)[a_iface];
      jweight = m_1inHomDirichWeight[a_dix](  jvof, 0)[  jface];

      idrop = m_0dirichletDropOrder[a_dix](a_ivof, 0)[a_iface];
      jdrop = m_1dirichletDropOrder[a_dix](  jvof, 0)[  jface];
    }
  else
    {
      jvof = m_1boundaryVofs[a_dix](a_ivof,0)[a_iface];
      jface = m_1boundaryFaceIndices[a_dix](a_ivof,0)[a_iface];

      a_istencil += m_1stencils[a_dix](a_ivof, 0)[a_iface];
      a_jstencil += m_0stencils[a_dix](  jvof, 0)[  jface];

      iweight = m_1inHomDirichWeight[a_dix](a_ivof, 0)[a_iface];
      jweight = m_0inHomDirichWeight[a_dix](  jvof, 0)[  jface];

      idrop = m_1dirichletDropOrder[a_dix](a_ivof, 0)[a_iface];
      jdrop = m_0dirichletDropOrder[a_dix](  jvof, 0)[  jface];
    }

  int jphase = 1-a_iphase;
  Real ibeta = a_beta[a_iphase];
  Real jbeta = a_beta[  jphase];
  Real factor = ibeta*iweight + jbeta*jweight;
  if (factor == 0)
    { // degenerate stencils, so have to give up
      pout() << "InterfaceJump::getFaceStencil - Degenerate stencils at "
             << a_ivof << endl;
      MayDay::Warning("Error in stencil to matrix conversion; matrix will be incorrect.");
    }
  else if (a_istencil.size()>0 && a_jstencil.size()>0)
    {
      a_istencil *=  jbeta*jweight/factor;
      a_jstencil *= -jbeta*iweight/factor;
    }
  else if (a_istencil.size()>0)
    {
      a_istencil.clear();
      a_jstencil.clear();
      EBISBox ebisBox;
      if (a_iphase == 0)
        {
          ebisBox = m_0ebisl[a_dix];
        }
      else
        {
          ebisBox = m_1ebisl[a_dix];
        }

      for (int idir=0; idir<CH_SPACEDIM; ++idir)
        {
          IntVect iv = a_ivof.gridIndex();

          iv[idir]+=1;

          int isign=1;
          int kdir = (idir+1) % CH_SPACEDIM;
          if (!m_domain.contains(iv) || ebisBox.isCovered(iv))
          {
            iv[idir]-=2;
            isign=-1;
          }
          if (!m_domain.contains(iv) || ebisBox.isCovered(iv))
          {
            iv[kdir]+=1;
          }
          if (!m_domain.contains(iv) || ebisBox.isCovered(iv))
          {
            iv[kdir]-=2;
          }
          if (!m_domain.contains(iv) || ebisBox.isCovered(iv))
          {
            iv[idir]+=2;
            isign=1;
          }
          if (!m_domain.contains(iv) || ebisBox.isCovered(iv))
          {
            iv[kdir]+=2;
          }
          if (!m_domain.contains(iv) || ebisBox.isCovered(iv))
          {
            iv[kdir]-=1;
            kdir=(idir+2) % CH_SPACEDIM;
            iv[kdir]+=1;
          }
          if (!m_domain.contains(iv) || ebisBox.isCovered(iv))
          {
            iv[kdir]-=2;
          }
          if (!m_domain.contains(iv) || ebisBox.isCovered(iv))
          {
            iv[idir]-=2;
            isign=-1;
          }

          CH_assert(!ebisBox.isCovered(iv));
          RealVect normal = ebisBox.normal(a_ivof, a_iface);
          VolIndex otherVof(iv, 0);
          Real otherWeight = normal[idir]*isign/m_vectDx[idir];
          a_istencil.add(otherVof, otherWeight);
          Real thisWeight = -normal[idir]*isign/m_vectDx[idir];
          a_istencil.add(a_ivof, thisWeight);
        }
      // switch sign to account for inward-pointing normal
      a_istencil *= -1.;
    }
  else if (a_jstencil.size()>0)
    {
      a_istencil.clear();
      a_jstencil.clear();
      EBISBox ebisBox;
      if (jphase == 0)
        {
          ebisBox = m_0ebisl[a_dix];
        }
      else
        {
          ebisBox = m_1ebisl[a_dix];
        }

      for (int idir=0; idir<CH_SPACEDIM; ++idir)
        {
          IntVect iv = jvof.gridIndex();

          iv[idir]+=1;

          int isign=1;
          int kdir = (idir+1) % CH_SPACEDIM;
          if (!m_domain.contains(iv) || ebisBox.isCovered(iv))
          {
            iv[idir]-=2;
            isign=-1;
          }
          if (!m_domain.contains(iv) || ebisBox.isCovered(iv))
          {
            iv[kdir]+=1;
          }
          if (!m_domain.contains(iv) || ebisBox.isCovered(iv))
          {
            iv[kdir]-=2;
          }
          if (!m_domain.contains(iv) || ebisBox.isCovered(iv))
          {
            iv[idir]+=2;
            isign=1;
          }
          if (!m_domain.contains(iv) || ebisBox.isCovered(iv))
          {
            iv[kdir]+=2;
          }
          if (!m_domain.contains(iv) || ebisBox.isCovered(iv))
          {
            iv[kdir]-=1;
            kdir=(idir+2) % CH_SPACEDIM;
            iv[kdir]+=1;
          }
          if (!m_domain.contains(iv) || ebisBox.isCovered(iv))
          {
            iv[kdir]-=2;
          }
          if (!m_domain.contains(iv) || ebisBox.isCovered(iv))
          {
            iv[idir]-=2;
            isign=-1;
          }

          CH_assert(!ebisBox.isCovered(iv));
          RealVect normal = ebisBox.normal(jvof, jface);
          VolIndex otherVof(iv, 0);
          Real otherWeight = normal[idir]*isign/m_vectDx[idir];
          a_jstencil.add(otherVof, otherWeight);
          Real thisWeight = -normal[idir]*isign/m_vectDx[idir];
          a_jstencil.add(jvof, thisWeight);
        }
      a_jstencil *= (jbeta/ibeta);
    }
  else
    {
      MayDay::Warning("InterfaceJump::getFaceStencil - Empty stencil!");
    }
}

void InterfaceJump::computeScalarBV(const EBCellFAB       & a_P0Fab,
                                    const EBCellFAB       & a_P1Fab,
                                    const Vector<Real>    & a_beta,
                                    const DataIndex       & a_dataInd,
                                    BaseIVFAB<Real>       & a_pBv0Fab,
                                    BaseIVFAB<Real>       & a_pBv1Fab,
                                    BaseIVFAB<Real>       & a_dpDn0Fab,
                                    BaseIVFAB<Real>       & a_dpDn1Fab,
                                    bool                    a_homogeneous)
{
  CH_TIME("InterfaceJump::computeScalarBV");
  CH_assert(!m_isVectorJump);

  int nPhase = 2;

  // need EBISBox for boundary areas
  Vector<EBISBox> ebisBoxV;
  ebisBoxV.resize(nPhase);
  ebisBoxV[0] = m_0ebisl[a_dataInd];
  ebisBoxV[1] = m_1ebisl[a_dataInd];

  Vector<VoFIterator> vofItV;
  vofItV.resize(nPhase);
  for (int iphase=0; iphase<nPhase; iphase++)
    {
      Box         region    = m_dbl.get(a_dataInd);
      IntVectSet  ivsBndry  = ebisBoxV[iphase].boundaryIVS(region);
      EBGraph     ebGraph   = ebisBoxV[iphase].getEBGraph();
      vofItV[iphase].define(ivsBndry, ebGraph);
    }

  //phi
  Vector<const EBCellFAB*> pFabs(nPhase);
  pFabs[0] = &a_P0Fab;
  pFabs[1] = &a_P1Fab;

  //boundary phi
  Vector< BaseIVFAB<Real>* > pBvFabs(nPhase);
  pBvFabs[0] = &a_pBv0Fab;
  pBvFabs[1] = &a_pBv1Fab;

  //boundary dPhi/dN
  Vector< BaseIVFAB<Real>* > dpDnFabs(nPhase);
  dpDnFabs[0] = &a_dpDn0Fab;
  dpDnFabs[1] = &a_dpDn1Fab;

  //stencil
  Vector< BaseIVFAB< Vector<VoFStencil> >* > stencils(nPhase);
  stencils[0] = &(m_0stencils[a_dataInd]);
  stencils[1] = &(m_1stencils[a_dataInd]);

  //inhomogeneous weights
  Vector< BaseIVFAB< Vector<Real> >* > inHomDirichWeightFab(nPhase);
  inHomDirichWeightFab[0] = &(m_0inHomDirichWeight[a_dataInd]);
  inHomDirichWeightFab[1] = &(m_1inHomDirichWeight[a_dataInd]);

  //drop order
  Vector< BaseIVFAB< Vector<int> >* > dirichletDropOrderFab(nPhase);
  dirichletDropOrderFab[0] = &(m_0dirichletDropOrder[a_dataInd]);
  dirichletDropOrderFab[1] = &(m_1dirichletDropOrder[a_dataInd]);

  //interface jump values
  BaseIVFAB< Vector<Real> >& gDFab = m_scalarGD[a_dataInd];
  BaseIVFAB< Vector<Real> >& gNFab = m_scalarGN[a_dataInd];

  // vofs and face indexes for other fluid
  Vector< BaseIVFAB< Vector<VolIndex> >* > boundaryVofFabs(nPhase);
  boundaryVofFabs[0] = &(m_0boundaryVofs[a_dataInd]);
  boundaryVofFabs[1] = &(m_1boundaryVofs[a_dataInd]);

  Vector< BaseIVFAB< Vector<int> >* > boundaryFaceIndexFabs(nPhase);
  boundaryFaceIndexFabs[0] = &(m_0boundaryFaceIndices[a_dataInd]);
  boundaryFaceIndexFabs[1] = &(m_1boundaryFaceIndices[a_dataInd]);

  Real pBj;
  Vector<Real> pBi;
  Real dPdNj;
  Vector<Real> dPdNi;

  for (int iphase=0; iphase<nPhase; iphase++)
    {
      // index of the other phase; assume two phases
      int jphase = 1-iphase;

      for (vofItV[iphase].reset(); vofItV[iphase].ok();
           ++vofItV[iphase])
        {
          VolIndex ivof = vofItV[iphase]();
          int niFace = ebisBoxV[iphase].numFacePhase(ivof);
          pBi.resize(niFace);
          dPdNi.resize(niFace);

          for (int iface=0; iface<niFace; iface++)
            {
              VolIndex jvof = (*boundaryVofFabs[iphase])(ivof,0)[iface];
              int jface = (*boundaryFaceIndexFabs[iphase])(ivof,0)[iface];

              Real gD, gN;
              if (!a_homogeneous)
                {
                  if (iphase == 0)
                    { // boundary jumps are defined in phase 0
                      gN = gNFab(ivof,0)[iface];
                      gD = gDFab(ivof,0)[iface];
                    }
                  else
                    {
                      gN = gNFab(jvof,0)[jface];
                      gD = gDFab(jvof,0)[jface];
                    }
                }
              else
                {
                  gN = 0.0;
                  gD = 0.0;
                }

              int dropOrderi = (*dirichletDropOrderFab[iphase])(ivof,0)[iface];
              int dropOrderj = (*dirichletDropOrderFab[jphase])(jvof,0)[jface];
              Real zi = (*inHomDirichWeightFab[iphase])(ivof, 0)[iface];
              Real zj = (*inHomDirichWeightFab[jphase])(jvof, 0)[jface];
              if ( (Abs(a_beta[iphase]*zi + a_beta[jphase]*zj) < 1E-6) &&
                  (dropOrderi == -1 && dropOrderj == -1) )
                { //yargh, both sides of interface formed a stencil, BUT, they
                  // ended up using the same grid points, so crap on one of them
                  dropOrderi = 1;
                }

              if (dropOrderi == -1 && dropOrderj == -1)
                {
                  Real rN, rD;
                  Real qi, qj;

                  rN = 0.0;
                  rD = 0.0;
                  qi = 0.0;
                  qj = 0.0;
                  pBi[iface] = 0.0;
                  pBj = 0.0;
                  dPdNi[iface] = 0.0;
                  dPdNj = 0.0;

                  const VoFStencil& vsi = (*stencils[iphase])(ivof,0)[iface];
                  int stenSizei = vsi.size();
                  const VoFStencil& vsj = (*stencils[jphase])(jvof,0)[jface];
                  int stenSizej = vsj.size();

                  //phase i q's
                  for (int kthTerm = 0; kthTerm < stenSizei; kthTerm++)
                    {
                      const VolIndex& kthVof = vsi.vof(kthTerm);
                      Real weight = vsi.weight(kthTerm);
                      qi += weight*(*pFabs[iphase])(kthVof, 0, 0);
                    }

                  //phase j q's
                  for (int kthTerm = 0; kthTerm < stenSizej; kthTerm++)
                    {
                      const VolIndex& kthVof = vsj.vof(kthTerm);
                      Real weight = vsj.weight(kthTerm);
                      qj += weight*(*pFabs[jphase])(kthVof, 0, 0);
                    }

                  //calculate zeta^s
                  Real zetai = -zi;
                  Real betaZetai = zetai*a_beta[iphase];

                  Real zetaj = -zj;
                  Real betaZetaj = zetaj*a_beta[jphase];

                  //calculate r^s
                  rN = gN - a_beta[iphase]*qi - a_beta[jphase]*qj;
                  rD = gD;

                  //calculate p^{Bs}
                  if (iphase == 0)
                    {
                      invert2x2(betaZetai, betaZetaj, rN, rD,
                                pBi[iface], pBj);
                    }
                  else
                    {
                      invert2x2(betaZetaj, betaZetai, rN, rD,
                                pBj, pBi[iface]);
                    }

                  //calculate dPdn^s
                  normalDeriv(qi, zetai, pBi[iface], dPdNi[iface]);
                  normalDeriv(qj, zetaj, pBj       , dPdNj);

                  if (iphase == 0)
                    {
                      dPdNi[iface] *= -1.0;
                    }
                  else
                    {
                      dPdNj *= -1.0;
                    }
                }
              else if (dropOrderi == -1)
                { // compute phi and dPhidN at the boundary just using phase i
                  // in fact, just ignore phi at boundary, since it is not used
                  // by EBAMRPoissonOp::applyOp anyways.
                  RealVect normal = ebisBoxV[iphase].normal(ivof, iface);
                  boundarySlope(&(dPdNi[iface]), normal, ivof,
                                *(pFabs[iphase]), ebisBoxV[iphase]);
                  if (iphase == 0)
                    {
                      dPdNj = (a_beta[iphase]*dPdNi[iface] + gN)
                              /a_beta[jphase];
                    }
                  else
                    {
                      dPdNj = (a_beta[iphase]*dPdNi[iface] - gN)
                              /a_beta[jphase];
                    }
                }
              else if (dropOrderj == -1)
                {
                  RealVect normal = ebisBoxV[jphase].normal(jvof, jface);
                  boundarySlope(&dPdNj, normal, jvof,
                                (*pFabs[jphase]), ebisBoxV[jphase]);
                  if (iphase == 0)
                    {
                      dPdNi[iface] = (a_beta[jphase]*dPdNj - gN)
                                     /a_beta[iphase];
                    }
                  else
                    {
                      dPdNi[iface] = (a_beta[jphase]*dPdNj + gN)
                                     /a_beta[iphase];
                    }
                }
              else
                { //we dropped order all the way through least squares,
                  // so set the phase 0 jump to zero and phi at the
                  // interface to the value in phase 0.
                  MayDay::Warning("InterfaceJump::computeScalarBV: Too coarse - dropping interface jump to no-stencil guess");
                  if (iphase == 0)
                    {
                      pBi[iface]   = (*pFabs[iphase])(ivof, 0);
                      pBj          = pBi[iface] + gD;
                      dPdNi[iface] = 0.0;
                      dPdNj        = gN/a_beta[jphase];
                    }
                  else
                    {
                      pBi[iface]   = (*pFabs[jphase])(jvof, 0) + gD;
                      pBj          = pBi[iface] - gD;
                      dPdNi[iface] = gN/a_beta[iphase];
                      dPdNj        = 0.0;
                    }
                }
            }

          //set return values
          (*pBvFabs[iphase])(ivof, 0) = 0.0;
          (*dpDnFabs[iphase])(ivof, 0) = 0.0;
          Real avgBndryArea = ebisBoxV[iphase].bndryArea(ivof);
          if (avgBndryArea > 0.0)
            {
              for (int iface=0; iface<niFace; iface++)
                {
                  Real bndryArea = ebisBoxV[iphase].bndryArea(ivof, iface);
                  ( *pBvFabs[iphase])(ivof, 0) += bndryArea*pBi[iface];
                  (*dpDnFabs[iphase])(ivof, 0) += bndryArea*dPdNi[iface];
                }
              ( *pBvFabs[iphase])(ivof, 0) /= avgBndryArea;
              (*dpDnFabs[iphase])(ivof, 0) /= avgBndryArea;
            }

          // since dpDn's were calculated with phase 0's normal,
          // we switch the sign of dpDN in phase 1
          if (iphase == 1)
            {
              (*dpDnFabs[iphase])(ivof, 0) *= -1.0;
            }
        }
    }
}

void InterfaceJump::computeVectorBV(const EBCellFAB       & a_P0Fab,
                                    const EBCellFAB       & a_P1Fab,
                                    const Vector<Real>    & a_beta,
                                    const DataIndex       & a_dataInd,
                                    BaseIVFAB<Real>       & a_pBv0Fab,
                                    BaseIVFAB<Real>       & a_pBv1Fab,
                                    BaseIVFAB<Real>       & a_dpDn0Fab,
                                    BaseIVFAB<Real>       & a_dpDn1Fab,
                                    bool                    a_homogeneous)
{
  CH_TIME("InterfaceJump::computeVectorBV");
  CH_assert(m_isVectorJump);

  int nPhase = 2;

  // need EBISBox for boundary areas
  Vector<EBISBox> ebisBoxV;
  ebisBoxV.resize(nPhase);
  ebisBoxV[0] = m_0ebisl[a_dataInd];
  ebisBoxV[1] = m_1ebisl[a_dataInd];

  Vector<VoFIterator> vofItV;
  vofItV.resize(nPhase);
  for (int iphase=0; iphase<nPhase; iphase++)
    {
      Box         region    = m_dbl.get(a_dataInd);
      IntVectSet  ivsBndry  = ebisBoxV[iphase].boundaryIVS(region);
      EBGraph     ebGraph   = ebisBoxV[iphase].getEBGraph();
      vofItV[iphase].define(ivsBndry, ebGraph);
    }

  //phi
  Vector<const EBCellFAB*> pFabs(nPhase);
  pFabs[0] = &a_P0Fab;
  pFabs[1] = &a_P1Fab;

  //boundary phi
  Vector< BaseIVFAB<Real>* > pBvFabs(nPhase);
  pBvFabs[0] = &a_pBv0Fab;
  pBvFabs[1] = &a_pBv1Fab;

  //boundary dPhi/dN
  Vector< BaseIVFAB<Real>* > dpDnFabs(nPhase);
  dpDnFabs[0] = &a_dpDn0Fab;
  dpDnFabs[1] = &a_dpDn1Fab;

  //stencil
  Vector< BaseIVFAB< Vector<VoFStencil> >* > stencils(nPhase);
  stencils[0] = &(m_0stencils[a_dataInd]);
  stencils[1] = &(m_1stencils[a_dataInd]);

  //inhomogeneous weights
  Vector< BaseIVFAB< Vector<Real> >* > inHomDirichWeightFab(nPhase);
  inHomDirichWeightFab[0] = &(m_0inHomDirichWeight[a_dataInd]);
  inHomDirichWeightFab[1] = &(m_1inHomDirichWeight[a_dataInd]);

  //drop order
  Vector< BaseIVFAB< Vector<int> >* > dirichletDropOrderFab(nPhase);
  dirichletDropOrderFab[0] = &(m_0dirichletDropOrder[a_dataInd]);
  dirichletDropOrderFab[1] = &(m_1dirichletDropOrder[a_dataInd]);

  //interface jump values
  BaseIVFAB< Vector<RealVect> >& gDFab = m_vectorGD[a_dataInd];
  BaseIVFAB< Vector<RealVect> >& gNFab = m_vectorGN[a_dataInd];

  // vofs and face indexes for other fluid
  Vector< BaseIVFAB< Vector<VolIndex> >* > boundaryVofFabs(nPhase);
  boundaryVofFabs[0] = &(m_0boundaryVofs[a_dataInd]);
  boundaryVofFabs[1] = &(m_1boundaryVofs[a_dataInd]);

  Vector< BaseIVFAB< Vector<int> >* > boundaryFaceIndexFabs(nPhase);
  boundaryFaceIndexFabs[0] = &(m_0boundaryFaceIndices[a_dataInd]);
  boundaryFaceIndexFabs[1] = &(m_1boundaryFaceIndices[a_dataInd]);

  RealVect pBj;
  Vector<RealVect> pBi;
  RealVect dPdNj;
  Vector<RealVect> dPdNi;

  for (int iphase=0; iphase<nPhase; iphase++)
    {
      // index of the other phase; assume two phases
      int jphase = 1-iphase;

      for (vofItV[iphase].reset(); vofItV[iphase].ok();
           ++vofItV[iphase])
        {
          VolIndex ivof = vofItV[iphase]();
          int niFace = ebisBoxV[iphase].numFacePhase(ivof);
          pBi.resize(niFace);
          dPdNi.resize(niFace);

          for (int iface=0; iface<niFace; iface++)
            {
              VolIndex jvof = (*boundaryVofFabs[iphase])(ivof,0)[iface];
              int jface = (*boundaryFaceIndexFabs[iphase])(ivof,0)[iface];

              RealVect gN, gD;
              if (!a_homogeneous)
                {
                  if (iphase == 0)
                    { // boundary jumps are defined in phase 0
                      gN = gNFab(ivof,0)[iface];
                      gD = gDFab(ivof,0)[iface];
                    }
                  else
                    {
                      gN = gNFab(jvof,0)[jface];
                      gD = gDFab(jvof,0)[jface];
                    }
                }
              else
                {
                  gN = RealVect::Zero;
                  gD = RealVect::Zero;
                }

              int dropOrderi = (*dirichletDropOrderFab[iphase])(ivof,0)[iface];
              int dropOrderj = (*dirichletDropOrderFab[jphase])(jvof,0)[jface];
              Real zi = (*inHomDirichWeightFab[iphase])(ivof, 0)[iface];
              Real zj = (*inHomDirichWeightFab[jphase])(jvof, 0)[jface];
              if ( (Abs(a_beta[iphase]*zi + a_beta[jphase]*zj) < 1E-6) &&
                  (dropOrderi == -1 && dropOrderj == -1) )
                { //yargh, both sides of interface formed a stencil, BUT, they
                  // ended up using the same grid points, so crap on one of them
                  dropOrderi = 1;
                }

              if (dropOrderi == -1 && dropOrderj == -1)
                {
                  RealVect rN, rD;
                  RealVect qi, qj;

                  rN = RealVect::Zero;
                  rD = RealVect::Zero;
                  qi = RealVect::Zero;
                  qj = RealVect::Zero;
                  pBi[iface] = RealVect::Zero;
                  pBj = RealVect::Zero;
                  dPdNi[iface] = RealVect::Zero;
                  dPdNj = RealVect::Zero;

                  const VoFStencil& vsi = (*stencils[iphase])(ivof,0)[iface];
                  int stenSizei = vsi.size();
                  const VoFStencil& vsj = (*stencils[jphase])(jvof,0)[jface];
                  int stenSizej = vsj.size();

                  //phase i q's
                  for (int idir=0; idir<SpaceDim; idir++)
                    {
                      for (int kthTerm = 0; kthTerm < stenSizei; kthTerm++)
                        {
                          const VolIndex& kthVof = vsi.vof(kthTerm);
                          Real weight = vsi.weight(kthTerm);
                          qi[idir] += weight*(*pFabs[iphase])(kthVof, idir, -1);
                        }

                      //phase j q's
                      for (int kthTerm = 0; kthTerm < stenSizej; kthTerm++)
                        {
                          const VolIndex& kthVof = vsj.vof(kthTerm);
                          Real weight = vsj.weight(kthTerm);
                          qj[idir] += weight*(*pFabs[jphase])(kthVof, idir, -1);
                        }
                    }

                  //calculate zeta^s
                  Real zetai = -zi;
                  Real betaZetai = zetai*a_beta[iphase];

                  Real zetaj = -zj;
                  Real betaZetaj = zetaj*a_beta[jphase];

                  //calculate r^s
                  rN = gN - a_beta[iphase]*qi - a_beta[jphase]*qj;
                  rD = gD;

                  for (int idir=0; idir<SpaceDim; idir++)
                    {
                      //calculate p^{Bs}
                      if (iphase == 0)
                        {
                          invert2x2(betaZetai, betaZetaj, rN[idir], rD[idir],
                                    pBi[iface][idir], pBj[idir]);
                        }
                      else
                        {
                          invert2x2(betaZetaj, betaZetai, rN[idir], rD[idir],
                                    pBj[idir], pBi[iface][idir]);
                        }

                      //calculate dPdn^s
                      normalDeriv(qi[idir], zetai, pBi[iface][idir], dPdNi[iface][idir]);
                      normalDeriv(qj[idir], zetaj,        pBj[idir]       , dPdNj[idir]);
                    }

                  if (iphase == 0)
                    {
                      dPdNi[iface] *= -1.0;
                    }
                  else
                    {
                      dPdNj *= -1.0;
                    }
                }
              else if (dropOrderi == -1)
                { // compute phi and dPhidN at the boundary just using phase i
                  // in fact, just ignore phi at boundary, since it is not used
                  // by EBAMRPoissonOp::applyOp anyways.
                  RealVect normal = ebisBoxV[iphase].normal(ivof, iface);
                  for (int idir=0; idir<SpaceDim; idir++)
                    {
                      boundarySlope(&(dPdNi[iface][idir]), normal, ivof,
                                    *(pFabs[iphase]), ebisBoxV[iphase], idir);
                    }
                  if (iphase == 0)
                    {
                      dPdNj = dPdNi[iface] + gN;
                    }
                  else
                    {
                      dPdNj = dPdNi[iface] - gN;
                    }
                }
              else if (dropOrderj == -1)
                {
                  RealVect normal = ebisBoxV[jphase].normal(jvof, jface);
                  for (int idir=0; idir<SpaceDim; idir++)
                    {
                      boundarySlope(&(dPdNj[idir]), normal, jvof,
                                    (*pFabs[jphase]), ebisBoxV[jphase], idir);
                    }
                  if (iphase == 0)
                    {
                      dPdNi[iface] = dPdNj - gN;
                    }
                  else
                    {
                      dPdNi[iface] = dPdNj + gN;
                    }
                }
              else
                { //we dropped order all the way through least squares,
                  // so set everything to zero
                  MayDay::Warning("InterfaceJump::computeScalarBV: Too coarse - dropping interface jump to no-stencil guess");
                  if (iphase == 0)
                    {
                      pBi[iface]   = RealVect::Zero;
                      pBj          = gD;
                      dPdNi[iface] = RealVect::Zero;
                      dPdNj        = gN;
                    }
                  else
                    {
                      pBi[iface]   = gD;
                      pBj          = RealVect::Zero;
                      dPdNi[iface] = gN;
                      dPdNj        = RealVect::Zero;
                    }
                }
            }

          //set return values
          for (int idir=0; idir<SpaceDim; idir++)
            {
              (*pBvFabs[iphase])(ivof, idir) = 0.0;
            }
          for (int iface=0; iface<niFace; iface++)
            {
              for (int idir=0; idir<SpaceDim; idir++)
                {
                  (*dpDnFabs[iphase])(ivof, idir) = 0.0;
                }
            }
          Real avgBndryArea = ebisBoxV[iphase].bndryArea(ivof);
          if (avgBndryArea > 0.0)
            {
              for (int idir=0; idir<SpaceDim; idir++)
                {
                  for (int iface=0; iface<niFace; iface++)
                    {
                      Real bndryArea = ebisBoxV[iphase].bndryArea(ivof, iface);
                      ( *pBvFabs[iphase])(ivof, idir) += bndryArea*pBi[iface][idir];
                      (*dpDnFabs[iphase])(ivof, idir) += bndryArea*dPdNi[iface][idir];
                    }
                  ( *pBvFabs[iphase])(ivof, 0) /= avgBndryArea;
                  (*dpDnFabs[iphase])(ivof, 0) /= avgBndryArea;
                }
            }

          // since dpDn's were calculated with phase 0's normal,
          // we switch the sign of dpDN in phase 1
          if (iphase == 1)
            {
              for (int idir=0; idir<SpaceDim; idir++)
                {
                  (*dpDnFabs[iphase])(ivof, idir) *= -1.0;
                }
            }
        }
    }
}

void InterfaceJump::invert2x2(const Real& a_betaZeta0,
                              const Real& a_betaZeta1,
                              const Real& a_rN,
                              const Real& a_rD,
                              Real&       a_pB0,
                              Real&       a_pB1)
{
  Real det = -(a_betaZeta1 + a_betaZeta0);
  a_pB1 = (1/det)*(-a_rN - (a_betaZeta0*a_rD));
  a_pB0 = (1/det)*(-a_rN + (a_betaZeta1*a_rD));

}

void  InterfaceJump::normalDeriv(const Real& a_qs,
                                 const Real& a_zetas,
                                 const Real& a_pBs,
                                 Real&       a_dPdns)
{
  a_dPdns = a_zetas*a_pBs + a_qs;
}

void InterfaceJump::findTangentPlane(RealVect&      a_t1,
                                     RealVect&      a_t2,
                                     const RealVect a_normal)
{

  Real modulus;

  Real c1[3];
  Real c2[3];
  Real cp[3];

  int  maxDir;
  Real maxVal = 0.0;

  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      if (fabs(a_normal[idir]) > maxVal)
        {
          maxDir = idir;
          maxVal = a_normal[idir];
        }

      c1[idir] = a_normal[idir];
      c2[idir] = a_normal[idir];
    }
  int otherDir;
  if (maxDir == SpaceDim - 1)
    {
      otherDir = 0;
    }
  else
    {
      otherDir = maxDir + 1;
    }
  c2[otherDir] += 1.0;

#if CH_SPACEDIM == 2
  {
    c1[2] = 0.0;
    c2[2] = 1.0;
  }
#endif
  crossProduct(cp,c1,c2);
  for (int idir = 0; idir < SpaceDim; ++ idir)
    {
      a_t1[idir] = cp[idir];
    }
  modulus =  euclidianNorm(a_t1);
  a_t1 *= 1/modulus;
  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      c1[idir] = a_normal[idir];
      c2[idir] = a_t1[idir];
    }
  if (SpaceDim == 2)
    {
      c1[2] = 0.0;
      c2[2] = 1.0;
    }
  crossProduct(cp,c1,c2);
  for (int idir = 0; idir < SpaceDim; ++ idir)
    {
      a_t2[idir] = cp[idir];
    }

  modulus =  euclidianNorm(a_t2);
  a_t2 *= 1/modulus;
}

void InterfaceJump::crossProduct(Real       a_cp[3],
                                 const Real a_c1[3],
                                 const Real a_c2[3])
{
  a_cp[0] =   a_c1[1]*a_c2[2] - a_c1[2]*a_c2[1];
  a_cp[1] = -(a_c1[0]*a_c2[2] - a_c1[2]*a_c2[0]);
  a_cp[2] =   a_c1[0]*a_c2[1] - a_c1[1]*a_c2[0];
}

Real InterfaceJump::euclidianNorm(const RealVect a_vect)
{
  Real retval = 0.0;
  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      retval += a_vect[idir]*a_vect[idir];
    }
  retval = sqrt(retval);
  return retval;
}
Real InterfaceJump::innerProduct(const RealVect& a_v1,
                                 const RealVect& a_v2)
{
  Real retval = 0.0;
  for (int idir = 0; idir < SpaceDim; ++idir)
    {
      retval += a_v1[idir]*a_v2[idir];
    }
  return retval;
}

void InterfaceJump::AInvTimesVelocity3D(RealVect&            a_answer,
                                        const RealVect&      a_velocity,
                                        const RealVect&      a_normal,
                                        const RealVect&      a_t1,
                                        const RealVect&      a_t2,
                                        const Vector<Real>&  a_nu,
                                        const Vector<Real>&  a_zeta)
{
  Real tangentFactor = 1.0/(a_nu[0]*a_zeta[0] +a_nu[1]*a_zeta[1]);
  Real normalFactor = 1.0/(a_zeta[0] +a_zeta[1]);

  RealVect Row0;
  Row0[0] = normalFactor*a_normal[0];
  Row0[1] = tangentFactor*a_t1[0];
  Row0[2] = tangentFactor*a_t2[0];
  a_answer[0] = innerProduct(Row0,a_velocity);

  RealVect Row1;
  Row1[0] = normalFactor*a_normal[1];
  Row1[1] = tangentFactor*a_t1[1];
  Row1[2] = tangentFactor*a_t2[1];
  a_answer[1] = innerProduct(Row1,a_velocity);

  RealVect Row2;
  Row2[0] = normalFactor*a_normal[2];
  Row2[1] = tangentFactor*a_t1[2];
  Row2[2] = tangentFactor*a_t2[2];
  a_answer[2] = innerProduct(Row2,a_velocity);

}

void InterfaceJump::AInvTimesVelocity2D(RealVect&            a_answer,
                                        const RealVect&      a_velocity,
                                        const RealVect&      a_normal,
                                        const RealVect&      a_t1,
                                        const Vector<Real>&  a_nu,
                                        const Vector<Real>&  a_zeta)
{
  Real tangentFactor = 1.0/(a_nu[0]*a_zeta[0] +a_nu[1]*a_zeta[1]);
  Real normalFactor = 1.0/(a_zeta[0] +a_zeta[1]);

  RealVect Row0;
  Row0[0] = normalFactor*a_normal[0];
  Row0[1] = tangentFactor*a_t1[0];
  RealVect Row1;
  Row1[0] = normalFactor*a_normal[1];
  Row1[1] = tangentFactor*a_t1[1];

  a_answer[0] = innerProduct(Row0,a_velocity);
  a_answer[1] = innerProduct(Row1,a_velocity);

}
void InterfaceJump::ATimesQ3D(RealVect&       a_answer,
                              const RealVect& a_q,
                              const RealVect& a_normal,
                              const RealVect& a_t1,
                              const RealVect& a_t2,
                              const Real      a_nu)
{
  a_answer[0] = innerProduct(a_normal,a_q);

  a_answer[1] = innerProduct(a_t1,a_q);
  a_answer[2] = innerProduct(a_t2,a_q);

  a_answer[1] *= a_nu;
  a_answer[2] *= a_nu;
}

void InterfaceJump::ATimesQ2D(RealVect&       a_answer,
                              const RealVect& a_q,
                              const RealVect& a_normal,
                              const RealVect& a_t1,
                              const Real      a_nu)
{
  a_answer[0] = innerProduct(a_normal,a_q);
  a_answer[1] = innerProduct(a_t1,a_q);

  a_answer[1] *= a_nu;
}
void InterfaceJump::ATimesQ(RealVect&       a_answer,
                            const RealVect& a_q,
                            const RealVect& a_normal,
                            const RealVect& a_t1,
                            const RealVect& a_t2,
                            const Real      a_nu)
{
#if CH_SPACEDIM == 3
  {
    ATimesQ3D(a_answer,a_q,a_normal,a_t1,a_t2,a_nu);
  }
#  elif CH_SPACEDIM == 2
  {
    ATimesQ2D(a_answer,a_q,a_normal,a_t1,a_nu);
  }
#endif
}
void InterfaceJump::AInvTimesVelocity(RealVect&            a_answer,
                                      const RealVect&      a_velocity,
                                      const RealVect&      a_normal,
                                      const RealVect&      a_t1,
                                      const RealVect&      a_t2,
                                      const Vector<Real>&  a_nu,
                                      const Vector<Real>&  a_zeta)
{
# if CH_SPACEDIM == 3
  AInvTimesVelocity3D(a_answer,
                      a_velocity,
                      a_normal,
                      a_t1,
                      a_t2,
                      a_nu,
                      a_zeta);
#elif CH_SPACEDIM == 2
  AInvTimesVelocity2D(a_answer,
                      a_velocity,
                      a_normal,
                      a_t1,
                      a_nu,
                      a_zeta);
#endif
}
#include "NamespaceFooter.H"

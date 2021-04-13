#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MFPoissonOp.H"
#include "BaseIVFactory.H"
#include "MFAliasFactory.H"
#include "MFLevelDataOps.H"
#include "CH_Timer.H"
#include <iomanip>
#include "NamespaceHeader.H"

class MFBC : public BaseEBBC
{
public:
  virtual LayoutData<BaseIVFAB<VoFStencil> >* getFluxStencil(int a_var)
  {
    return NULL;
  }
  virtual void define(const LayoutData<IntVectSet>& a_cfivs,
                      const Real&                   a_factor)
  {
  }

  MFBC()
  {
  }

  virtual ~MFBC()
  {
  }

  virtual void getEBFlux(Real&                       a_flux,
                         const VolIndex&             a_vof,
                         const LevelData<EBCellFAB>& a_phi,
                         const LayoutData<IntVectSet>& a_cfivs,
                         const DataIndex&            a_dit,
                         const RealVect&             a_probLo,
                         const RealVect&             a_dx,
                         const bool&                 a_useHomogeneous,
                         const Real&                 a_time,
                         const std::pair<int,Real>*  a_cacheHint=0 );
  virtual void applyEBFlux(EBCellFAB&                    a_lphi,
                           const EBCellFAB&              a_phi,
                           VoFIterator&                  a_vofit,
                           const LayoutData<IntVectSet>& a_cfivs,
                           const DataIndex&              a_dit,
                           const RealVect&               a_probLo,
                           const RealVect&               a_dx,
                           const Real&                   a_factor,
                           const bool&                   a_useHomogeneous,
                           const Real&                   a_time)
  {
  }
};

void MFBC::getEBFlux(Real&                       a_flux,
                     const VolIndex&             a_vof,
                     const LevelData<EBCellFAB>& a_phi,
                     const LayoutData<IntVectSet>& a_cfivs,
                     const DataIndex&            a_dit,
                     const RealVect&             a_probLo,
                     const RealVect&             a_dx,
                     const bool&                 a_useHomogeneous,
                     const Real&                 a_time,
                     const std::pair<int,Real>*  a_cacheHint )
{
  CH_assert(a_useHomogeneous); // this BC should only be called for homogeneous solving.
  a_flux = 0; // pretty simple, homogeneous flux at boundary is zero.

}

MFPoissonOp::~MFPoissonOp()
{
  for (int i=0; i<m_alias.size(); i++)
    {
      delete m_alias[i];
    }
  for (int i=0; i<m_ebops.size(); i++)
    {
      delete m_ebops[i];
    }
}
/**
   \name LinearOp functions */
/*@{*/


/** full define function for AMRLevelOp with both coarser and finer levels */
void MFPoissonOp::define(const MFIndexSpace&                         a_mfis,
                         int                                         a_ncomp,
                         const DisjointBoxLayout&                    a_grids,
                         const DisjointBoxLayout&                    a_gridsCoarMG,
                         const bool&                                 a_hasMGObjects,
                         const bool&                                 a_layoutChanged,
                         const DisjointBoxLayout&                    a_gridsFiner,
                         const DisjointBoxLayout&                    a_gridsCoarser,
                         const RealVect&                             a_dxLevel,
                         int                                         a_refRatio,
                         int                                         a_refRatioFiner,
                         const ProblemDomain&                        a_domain,
                         const Vector<RefCountedPtr<BaseDomainBC> >& a_bc,
                         const IntVect&                              a_ghostPhi,
                         const IntVect&                              a_ghostRHS,
                         bool                                        hasCoarser,
                         bool                                        hasFiner,
                         const Vector<Real>&                         a_alpha,
                         const Vector<Real>&                         a_beta)
{
  CH_TIME("MFPoissonOp::define");
  m_phases = a_mfis.numPhases();

  m_ghostPhi = a_ghostPhi;
  m_ghostRHS = a_ghostRHS;
  m_alpha = a_alpha;
  m_aCoef = a_alpha;
  m_beta  = a_beta;
  m_bCoef = a_beta;

  m_ncomp = a_ncomp;
  CH_assert(m_ncomp == 1 || m_ncomp == CH_SPACEDIM);
  m_refToCoarser = a_refRatio;
  m_ebops.resize(m_phases);

  m_alias.resize(4);
  int nghost = 4;
  m_dx = a_dxLevel;
  RefCountedPtr<BaseEBBC> mfbc( new MFBC );
  ProblemDomain cDomain = a_domain, fDomain = a_domain;
  RealVect dxCoarse = a_dxLevel * a_refRatio;
  RealVect dxFine   = a_dxLevel / a_refRatioFiner;

  m_domain = a_domain;
  cDomain.coarsen(a_refRatio);
  fDomain.refine(a_refRatioFiner);
  m_alias[0] = new LevelData<EBCellFAB>;
  m_alias[1] = new LevelData<EBCellFAB>;
  m_alias[2] = new LevelData<EBCellFAB>;
  m_alias[3] = new LevelData<EBCellFAB>;

  Vector<EBISLayout> layouts(m_phases);
  Vector<int> ncomp(m_phases, m_ncomp);
  m_relax = 0;
  RealVect origin = a_mfis.EBIS(0)->getOrigin();
  bool isVectorJump = (m_ncomp == CH_SPACEDIM);
  m_jump.define(a_mfis, a_domain, a_grids, a_dxLevel, origin, isVectorJump);

  //for now set the ghost cells to a default.  need
  //to get this into the interface of the operator.
  //MayDay::Warning("multigrid coarsened versions of ebislayouts and dbls not defined");
  DisjointBoxLayout dblCoarMG = a_gridsCoarMG;

  ProblemDomain domainCoarMG = coarsen(a_domain, 2);
  bool layoutChanged = a_layoutChanged;
  bool hasMGObjects =  a_hasMGObjects;


  for (int i=0; i<m_phases; i++)
    {
      const EBIndexSpace* ebisPtr = a_mfis.EBIS(i);
      EBLevelGrid eblg, eblgFine, eblgCoar;
      eblg = EBLevelGrid(a_grids,a_domain, nghost, ebisPtr);

      RealVect dxCoar = a_dxLevel*a_refRatio;  //not technically correct

      BaseIVFactory<Real> factory(eblg.getEBISL());
      m_boundaryN[i].define(a_grids, m_ncomp, IntVect::Zero, factory);
      m_boundaryD[i].define(a_grids, m_ncomp, IntVect::Zero, factory);

      RefCountedPtr<EBQuadCFInterp> quadCFI;
      if (hasFiner)
        {
          eblgFine = EBLevelGrid(a_gridsFiner, fDomain, nghost, ebisPtr);
        }
      if (hasCoarser)
        {
          eblgCoar = EBLevelGrid(a_gridsCoarser, cDomain, nghost, ebisPtr);
          quadCFI = RefCountedPtr<EBQuadCFInterp>
            (new EBQuadCFInterp(a_grids, a_gridsCoarser,
                                eblg.getEBISL(), eblgCoar.getEBISL(),
                                cDomain, a_refRatio, 1,
                                *eblg.getCFIVS(), ebisPtr));
        }
      //layoutChanged = true;
      //hasMGObjects = false; //hack.   to get things running again
      EBLevelGrid eblgCoarMG;
      if (hasMGObjects)
        {
          eblgCoarMG = EBLevelGrid(dblCoarMG, domainCoarMG, nghost, a_mfis.EBIS(i));
        }
      layouts[i] = eblg.getEBISL();
      m_ebops[i] = new EBAMRPoissonOp(eblgFine, eblg, eblgCoar, eblgCoarMG,
                                      quadCFI, a_bc[i], mfbc, a_dxLevel, dxCoar,
                                      origin,a_refRatioFiner,a_refRatio,
                                      hasFiner, hasCoarser, hasMGObjects, layoutChanged,
                                      40, 0, m_alpha[i], m_beta[i], m_ghostPhi, m_ghostRHS);


    }
  MFCellFactory* factory = new MFCellFactory(layouts, ncomp);
  m_tmp.define( a_grids, m_ncomp, m_ghostRHS, *factory);
  m_weights.define( a_grids, m_ncomp, m_ghostPhi, *factory);
  RefCountedPtr<DataFactory<MFCellFAB> > fac(factory);

  m_ops.define(fac);// fac gets deleted by RefCountedPtr.

  m_colors.resize(EBAMRPO_NUMSTEN);
#if CH_SPACEDIM==2
  m_colors[0] = IntVect::Zero;//(0,0)
  m_colors[1] = IntVect::Unit;//(1,1)
  m_colors[2] = IntVect::Zero + BASISV(1);//(0,1)
  m_colors[3] = IntVect::Zero + BASISV(0);//(1,0)
#elif CH_SPACEDIM==3
  m_colors[0] = IntVect::Zero;//(0,0,0)
  m_colors[1] = IntVect::Zero + BASISV(0) + BASISV(1);//(1,1,0)
  m_colors[2] = IntVect::Zero + BASISV(1) + BASISV(2);//(0,1,1)
  m_colors[3] = IntVect::Zero + BASISV(0) + BASISV(2);//(1,0,1)
  m_colors[4] = IntVect::Zero + BASISV(1);//(0,1,0)
  m_colors[5] = IntVect::Zero + BASISV(0);//(1,0,0)
  m_colors[6] = IntVect::Zero + BASISV(2);//(0,0,1)
  m_colors[7] = IntVect::Unit;//(1,1,1)
#endif

  /*
    Uncomment this in order to dump a matrix representation of the
    elliptic operator.
  */
  // dumpStencilMatrix();
}

void MFPoissonOp::setAlphaAndBeta(const Real& a_alpha,
                                  const Real& a_beta)
{
  CH_TIME("MFPoissonOp::setAlphaAndBeta");
  // CH_assert(a_alpha.size()==m_phases && a_beta.size()==m_phases);
  for (int iphase=0; iphase<m_phases; iphase++)
    {
      m_ebops[iphase]->setAlphaAndBeta(a_alpha, a_beta);
      m_alpha[iphase] = a_alpha*m_aCoef[iphase];
      m_beta[iphase]  = a_beta*m_bCoef[iphase];
    }
}

void MFPoissonOp::setTime(Real a_time)
{
  CH_TIME("MFPoissonOp::setTime");
  for (int iphase=0; iphase<m_phases; iphase++)
    {
      m_ebops[iphase]->setTime(a_time);
    }
  m_jump.resetJump(m_aCoef, m_bCoef, a_time);
}

void MFPoissonOp::diagonalScale(LevelData<MFCellFAB>& a_rhs)
{
  CH_TIME("MFPoissonOp::diagonalScale");
  MFLevelDataOps::kappaWeight(a_rhs);
}

void MFPoissonOp::applyOpNoBoundary(LevelData<MFCellFAB>&       a_opPhi,
                                    const LevelData<MFCellFAB>& a_phi)
{
  CH_TIME("MFPoissonOp::applyOpNoBoundary");
  for (int iphase=0; iphase<m_phases; iphase++)
    {
      aliasMF(*m_alias[0], iphase, a_opPhi);
      aliasMF(*m_alias[1], iphase, a_phi);
      m_ebops[iphase]->applyOpNoBoundary(*m_alias[0], *m_alias[1]);
    }
}

Real MFPoissonOp::totalBoundaryFlux(int                          a_phase,
                                    const LevelData<MFCellFAB>&  a_phi,
                                    Real                         a_factor,
                                    bool                         a_divideByTotalArea)
{
  CH_TIME("MFPoissonOp::totalBoundaryFlux");
  computeBoundaryN(a_phi, false);

  Real flux = 0.;
  Real area = 0.;
  EBLevelGrid eblg = m_ebops[a_phase]->getEBLG();
  EBISLayout ebisl = eblg.getEBISL();
  for (DataIterator dit = eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      BaseIVFAB<Real>& fluxIVFab = m_boundaryN[a_phase][dit];
      EBISBox ebisBox = ebisl[dit()];
      const Box& box = eblg.getDBL().get(dit());
      const EBGraph& ebgraph = ebisBox.getEBGraph();
      IntVectSet ivs = ebisBox.boundaryIVS(box);
      VoFIterator vofit(ivs, ebgraph);
      for (vofit.reset(); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real bndryArea = ebisBox.bndryArea(vof);
          area += bndryArea;
          flux += -bndryArea*fluxIVFab(vof, 0);
        }
    }
#ifdef CH_MPI
  {
    Real recv;
    int result;
    result = MPI_Allreduce(&flux, &recv, 1, MPI_CH_REAL,
                           MPI_SUM, Chombo_MPI::comm);
    flux = recv;
  }
#endif
  flux *= m_bCoef[a_phase]*a_factor;
  if (a_divideByTotalArea)
    {
      flux /= area;
    }
  return flux;
}

void MFPoissonOp::getBoundaryValues(LevelData<MFCellFAB>&  a_phi,
                                    LevelData<MFCellFAB>&  a_dPhi_dN,
                                    Real                   a_invalidVal)
{
  // first define a LevelData<FArrayBox> into which to copy
  // the values in a_dataPtr
  const DisjointBoxLayout& grids = m_boundaryD[0].getBoxes();
  int numComp = m_boundaryD[0].nComp();
  IntVect ghostVect = m_boundaryD[0].ghostVect();
  Vector<int> ncompVect(m_phases, numComp);
  Vector<EBISLayout> ebislVect(m_phases);
  for (int iphase = 0; iphase < m_phases; iphase++)
    {
      EBLevelGrid eblg = m_ebops[iphase]->getEBLG();
      ebislVect[iphase] = eblg.getEBISL();
    }
  MFCellFactory mfcfactory(ebislVect, ncompVect);
  a_phi.define(grids, numComp, ghostVect, mfcfactory);
  a_dPhi_dN.define(grids, numComp, ghostVect, mfcfactory);

  // now copy values in a_dataPtr into ldf
  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      for (int iphase = 0; iphase < m_phases; iphase++)
        {
          BaseIVFAB<Real>&     phiIV = m_boundaryD[iphase][dit()];
          BaseIVFAB<Real>& dPhi_dNIV = m_boundaryN[iphase][dit()];

          EBCellFAB&     phiEBCF =     a_phi[dit()].getPhase(iphase);
          EBCellFAB& dPhi_dNEBCF = a_dPhi_dN[dit()].getPhase(iphase);

          // set everything to a bogus value; all interface data
          // will be overwritten below
          phiEBCF.setVal(a_invalidVal);
          dPhi_dNEBCF.setVal(a_invalidVal);

          // now set values in EBCellFAB to values in boundary data.
          const EBGraph& ebgraph = phiIV.getEBGraph();
          const IntVectSet& ivsIrreg = phiIV.getIVS();
          for (VoFIterator vofit(ivsIrreg, ebgraph) ; vofit.ok(); ++vofit)
            {
              const VolIndex& vof = vofit();
              for (int icomp = 0; icomp < numComp; icomp++)
                {
                      phiEBCF(vof, icomp) =     phiIV(vof, icomp);
                  dPhi_dNEBCF(vof, icomp) = dPhi_dNIV(vof, icomp);
                }
            } // end loop over irregular vofs
        } // end loop over phases
    } // end loop over grids
}

Real MFPoissonOp::exactBoundaryFlux(int                         a_phase,
                                    RefCountedPtr<BaseBCValue>  a_flxVal,
                                    RealVect&                   a_origin,
                                    const Real&                 a_time)
{
  CH_TIME("MFPoissonOp::exactBoundaryFlux");
  Real flux = 0.;
  Real area = 0.;
  EBLevelGrid eblg = m_ebops[a_phase]->getEBLG();
  EBISLayout ebisl = eblg.getEBISL();
  for (DataIterator dit = eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      EBISBox ebisBox = ebisl[dit()];
      const Box& box = eblg.getDBL().get(dit());
      const EBGraph& ebgraph = ebisBox.getEBGraph();
      IntVectSet ivs = ebisBox.boundaryIVS(box);
      VoFIterator vofit(ivs, ebgraph);
      for (vofit.reset(); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          int nface = ebisBox.numFacePhase(vof);
          for (int iface=0; iface<nface; iface++)
            {
              RealVect normal = ebisBox.normal(vof, iface);
              RealVect loc = vof.gridIndex();
              loc += 0.5;
              loc += ebisBox.bndryCentroid(vof, iface);
              loc *= m_dx;
              loc += a_origin;

              Real bndryArea = ebisBox.bndryArea(vof, iface);
              area += bndryArea;
              flux += bndryArea*a_flxVal->value(loc, normal, a_time, 0);
            }
        }
    }
  flux *= m_bCoef[a_phase]/area;
  return flux;
}

void MFPoissonOp::dumpStencilMatrix()
{
  map<StencilIndex, Real, StencilIndexComparator> mapper;

  // Get the Level and EB data
  const IntVect& domainSize = m_domain.size();

  pout() << "% Matlab output for VoFStencil:" << endl;
  pout() << "\tN = " << domainSize << ";" << endl;
  pout() << "\th = " << m_dx << ";" << endl;

  // loop through phases adding non-covered VoFs
  for (int iphase=0; iphase<m_phases; ++iphase)
    {
      const EBLevelGrid& eblg = m_ebops[iphase]->getEBLG();
      const DisjointBoxLayout& dbl = eblg.getDBL();
      for (DataIterator dit=dbl.dataIterator(); dit.ok(); ++dit)
        {
          Box box = dbl[dit()];
          IntVectSet ivsall(box);
          const EBISBox& ebisBox = eblg.getEBISL()[dit()];
          const EBGraph& ebGraph = ebisBox.getEBGraph();
          for (VoFIterator vofit(ivsall, ebGraph); vofit.ok(); ++vofit)
            {
              VolIndex srcVof(vofit());
              IntVect srciv = srcVof.gridIndex();
              if (!ebisBox.isCovered(srciv))
                {
                  StencilIndex si(srcVof, iphase);
                  mapper[si] = ebisBox.volFrac(srcVof);
                }
            }
        }
    }

  // all the VoFs have been added to the map, and it is properly sorted.
  //  so, we just iterate through it setting the map_values sequentially.
  char charstr[100];
  map<StencilIndex, Real>::iterator mit;
  int imap = 0;
  for (mit = mapper.begin(); mit != mapper.end(); ++mit)
    {
      Real kappa = mit->second;
      IntVect iv = mit->first.vof().gridIndex();
      mit->second = ++imap;
#if CH_SPACEDIM==2
      sprintf(charstr, "M_%03d(%03d,:) = [ %4d %4d %1d %22.16e ];",
              domainSize[0], (int) mit->second,
              iv[0], iv[1], mit->first.phase(),
              kappa);
#elif CH_SPACEDIM==3
      sprintf(charstr, "M_%03d(%03d,:) = [ %4d %4d %4d %1d %22.16e ];",
              domainSize[0], (int) mit->second,
              iv[0], iv[1], iv[2], mit->first.phase(),
              kappa);
#endif
      string str(charstr);
      pout() << str << endl;
    }

  // loop through phases outputting stencil weights
  int comp = 0;
  for (int iphase=0; iphase<m_phases; ++iphase)
    {
      const EBLevelGrid& eblg = m_ebops[iphase]->getEBLG();
      const DisjointBoxLayout& dbl = eblg.getDBL();
      for (DataIterator dit=dbl.dataIterator(); dit.ok(); ++dit)
        {
          Box box = dbl[dit()];
          IntVectSet ivsall(box);
          const EBISBox& ebisBox = eblg.getEBISL()[dit()];
          const EBGraph& ebGraph = ebisBox.getEBGraph();
          for (VoFIterator vofit(ivsall, ebGraph); vofit.ok(); ++vofit)
            {
              VolIndex srcVof(vofit());
              IntVect srciv = srcVof.gridIndex();
              if (!ebisBox.isCovered(srciv))
                {
                  StencilIndex si(srcVof, iphase);
                  map<StencilIndex, Real>::iterator srcit = mapper.find(si);
                  if (srcit == mapper.end())
                    {
                      pout() << "could not find " << srcVof << ", phase=" << iphase << endl;
                      MayDay::Error("Source VoF not found in map");
                    }
                  else
                    {
                      // mapping for srcVof
                      int is = (int) srcit->second;
                      // face stencil
                      VoFStencil isten = VoFStencil();
                      m_ebops[iphase]->getDivFStencil(isten, srcVof, dit(), true);
                      {
                        VoFStencil domainFluxStencil = VoFStencil();
                        m_ebops[iphase]->getDomainFluxStencil(domainFluxStencil,
                                                              srcVof,
                                                              comp,
                                                              dit());
                        isten += domainFluxStencil;
                      }
                      VoFStencil jsten = VoFStencil();
                      if (ebisBox.isIrregular(srciv))
                        {
                          for (int iface=0; iface<ebisBox.numFacePhase(srcVof); ++iface)
                            {
                              Real boundaryArea = ebisBox.bndryArea(srcVof, iface);
                              if (boundaryArea > 0)
                                { // needed since interface jump only has data for
                                  //  irregular VoFs with non-zero boundary areas
                                  // boundary stencil per face
                                  VoFStencil iisten, jjsten;
                                  m_jump.getFaceStencil(iisten, jjsten, m_beta,
                                                        srcVof, dit(), iphase,
                                                        iface);
                                  // properly scale irregular stencils
                                  Real scale = boundaryArea/m_dx[0];
                                  iisten *= scale;
                                  jjsten *= scale;
                                  // add EB stencil to face stencil
                                  isten += iisten;
                                  jsten += jjsten;
                                }
                            }
                        }
                      // scale stencil by beta
                      isten *= m_beta[iphase];
                      // add alpha. the alpha weight is hard-wired to be kappa here,
                      //  since I do not have access to the EB's alpha weighting.
                      if (m_alpha[iphase] != 0)
                        {
                          isten.add(srcVof,
                                    m_alpha[iphase]*ebisBox.volFrac(srcVof),
                                    comp);
                        }
                      for (int i=0; i<isten.size(); ++i)
                        {
                          VolIndex dstVof(isten.vof(i));
                          StencilIndex di(dstVof, iphase);
                          map<StencilIndex, Real>::iterator dstit = mapper.find(di);
                          if (dstit == mapper.end())
                            {
                              IntVectSet ivsi = ebisBox.boundaryIVS(box);
                              if (ivsi.contains(dstVof.gridIndex()))
                                {
                                  pout() << "could not find " << dstVof << ", phase=" << iphase << endl;
                                  MayDay::Error("Destination VoF not found in map");
                                }
                            }
                          else if (isten.weight(i) != 0.)
                            {
                              // mapping for dstVof
                              int id = (int) dstit->second;
                              sprintf(charstr,"L_%03d(%03d, %03d) = %22.16e;",
                                      domainSize[0], is, id, isten.weight(i));
                              string str(charstr);
                              pout() << str << endl;
                            }
                        }
                      // add stencil from other phase
                      // scale stencil by beta
                      jsten *= m_beta[iphase];
                      int jphase = 1-iphase;
                      for (int i=0; i<jsten.size(); ++i)
                        {
                          VolIndex dstVof(jsten.vof(i));
                          StencilIndex di(dstVof, jphase);
                          map<StencilIndex, Real>::iterator dstit = mapper.find(di);
                          if (dstit == mapper.end())
                            {
                              const EBLevelGrid& eblgj = m_ebops[jphase]->getEBLG();
                              const EBISBox& ebisBoxj = eblgj.getEBISL()[dit()];
                              IntVectSet ivsj = ebisBoxj.boundaryIVS(box);
                              if (ivsj.contains(dstVof.gridIndex()))
                                {
                                  pout() << "could not find " << dstVof << ", phase=" << jphase << endl;
                                  MayDay::Error("Destination VoF not found in map");
                                }
                            }
                          else if (jsten.weight(i) != 0.)
                            {
                              // mapping for dstVof
                              int id = (int) dstit->second;
                              sprintf(charstr,"L_%03d(%03d, %03d) = %22.16e;",
                                      domainSize[0], is, id, jsten.weight(i));
                              string str(charstr);
                              pout() << str << endl;
                            }
                        }
                    }
                }
            }
        }
    }
}

void MFPoissonOp::dumpReferenceStencilMatrix()
{
  map<StencilIndex, Real, StencilIndexComparator> mapper;

  // Get the Level and EB data
  const IntVect& domainSize = m_domain.size();

  pout() << "% Matlab output for VoFStencil:" << endl;
  pout() << "\tN = " << domainSize << ";" << endl;
  pout() << "\th = " << m_dx << ";" << endl;

  // loop through phases adding non-covered VoFs
  for (int iphase=0; iphase<m_phases; ++iphase)
    {
      const EBLevelGrid& eblg = m_ebops[iphase]->getEBLG();
      const DisjointBoxLayout& dbl = eblg.getDBL();
      for (DataIterator dit=dbl.dataIterator(); dit.ok(); ++dit)
        {
          Box box = dbl[dit()];
          IntVectSet ivsall(box);
          const EBISBox& ebisBox = eblg.getEBISL()[dit()];
          const EBGraph& ebGraph = ebisBox.getEBGraph();
          for (VoFIterator vofit(ivsall, ebGraph); vofit.ok(); ++vofit)
            {
              VolIndex srcVof(vofit());
              IntVect srciv = srcVof.gridIndex();
              if (!ebisBox.isCovered(srciv))
                {
                  StencilIndex si(srcVof, iphase);
                  mapper[si] = ebisBox.volFrac(srcVof);
                }
            }
        }
    }

  // all the VoFs have been added to the map, and it is properly sorted.
  //  so, we just iterate through it setting the map_values sequentially.
  char charstr[100];
  map<StencilIndex, Real>::iterator mit;
  int imap = 0;
  for (mit = mapper.begin(); mit != mapper.end(); ++mit)
    {
      Real kappa = mit->second;
      IntVect iv = mit->first.vof().gridIndex();
      mit->second = ++imap;
#if CH_SPACEDIM==2
      sprintf(charstr, "Mref_%03d(%03d,:) = [ %4d %4d %1d %22.16e ];",
              domainSize[0], (int) mit->second,
              iv[0], iv[1], mit->first.phase(),
              kappa);
#elif CH_SPACEDIM==3
      sprintf(charstr, "Mref_%03d(%03d,:) = [ %4d %4d %4d %1d %22.16e ];",
              domainSize[0], (int) mit->second,
              iv[0], iv[1], iv[2], mit->first.phase(),
              kappa);
#endif
      string str(charstr);
      pout() << str << endl;
    }
  // pout() << "N2 = " << imap << ";" << endl;;

  // loop through phases outputting stencil weights
  LevelData<MFCellFAB>& lhs = m_tmp;
  int nvar = 1;
  Vector<int> nComp(m_phases,nvar);
  Vector<EBISLayout> levEbislv;
  levEbislv.resize(m_phases);
  for (int iphase=0; iphase<m_phases; ++iphase)
    {
      levEbislv[iphase] = m_ebops[iphase]->getEBLG().getEBISL();
    }
  MFCellFactory factory(levEbislv, nComp);
  LevelData<MFCellFAB> phi(m_ebops[0]->getEBLG().getDBL(),
                           nvar, m_ghostPhi, factory);

  for (int iphase=0; iphase<m_phases; ++iphase)
    {
      const EBLevelGrid& eblg = m_ebops[iphase]->getEBLG();
      const DisjointBoxLayout& dbl = eblg.getDBL();
      for (DataIterator dit=dbl.dataIterator(); dit.ok(); ++dit)
        {
          Box box = dbl[dit()];
          IntVectSet ivsall(box);
          const EBISBox& ebisBox = eblg.getEBISL()[dit()];
          const EBGraph& ebGraph = ebisBox.getEBGraph();
          EBCellFAB& data = phi[dit()].getPhase(iphase);
          for (VoFIterator vofit(ivsall, ebGraph); vofit.ok(); ++vofit)
            {
              VolIndex srcVof(vofit());
              IntVect srciv = srcVof.gridIndex();
              if (!ebisBox.isCovered(srciv))
                {
                  StencilIndex si(srcVof, iphase);
                  map<StencilIndex, Real>::iterator srcit = mapper.find(si);
                  if (srcit == mapper.end())
                    {
                      pout() << "could not find " << srcVof << ", phase=" << iphase << endl;
                      MayDay::Error("Source VoF not found in map");
                    }
                  else
                    {
                      // mapping for srcVof
                      int is = (int) srcit->second;
                      // calculate weights for srcVof
                      setToZero(lhs);
                      setToZero(phi);
                      data(srcVof, 0) = 1;
                      applyOp(lhs, phi, true);
                      // loop through vofs looking for non-zero entries
                      for (int jphase=0; jphase<m_phases; ++jphase)
                        {
                          const EBLevelGrid& eblgL = m_ebops[jphase]->getEBLG();
                          const DisjointBoxLayout& dblL = eblgL.getDBL();
                          for (DataIterator ditL=dblL.dataIterator(); ditL.ok(); ++ditL)
                            {
                              Box boxL = dblL[ditL()];
                              IntVectSet ivsL(boxL);
                              const EBISBox& ebisBoxL = eblgL.getEBISL()[ditL()];
                              const EBGraph& ebGraphL = ebisBoxL.getEBGraph();
                              EBCellFAB& dataL = lhs[ditL()].getPhase(jphase);
                              for (VoFIterator vofitL(ivsL, ebGraphL); vofitL.ok(); ++vofitL)
                                {
                                  VolIndex lVof(vofitL());
                                  IntVect liv = lVof.gridIndex();
                                  if ( (!ebisBox.isCovered(liv)) &&
                                       (dataL(lVof, 0) != 0.) )
                                    {
                                      StencilIndex li(lVof, jphase);
                                      map<StencilIndex, Real>::iterator lit = mapper.find(li);
                                      if (lit == mapper.end())
                                        {
                                          pout() << "could not find " << lVof
                                                 << ", phase=" << jphase << endl;
                                          MayDay::Error("Destination VoF not found in map");
                                        }
                                      else
                                        {
                                          // mapping for lVof
                                          int il = (int) lit->second;
                                          sprintf(charstr,"Lref_%03d(%03d, %03d) = %22.16e;",
                                                  domainSize[0], il, is, dataL(lVof, 0));
                                          string str(charstr);
                                          pout() << str << endl;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void MFPoissonOp::setJump(const Real& a_gD, const Real& a_gN)
{
  m_jump.setJump(a_gD, a_gN);
}

void MFPoissonOp::setJump(const RealVect& a_gD, const RealVect& a_gN)
{
  CH_TIME("MFPoissonOp::setJump(const RealVect)");
  m_jump.setJump(a_gD, a_gN);
}

void MFPoissonOp::setJump(const Vector< RefCountedPtr<BaseBCValue> >& a_phiValVect,
                          const Vector< RefCountedPtr<BaseBCValue> >& a_flxValVect)
{
  CH_TIME("MFPoissonOp::setJump(BaseBCValue)");
  m_jump.setJump(a_phiValVect, a_flxValVect, m_aCoef, m_bCoef);
}

void MFPoissonOp::residual( LevelData<MFCellFAB>&        a_lhs,
                            const LevelData<MFCellFAB>&  a_phi,
                            const LevelData<MFCellFAB>&  a_rhs,
                            bool                         a_homogeneous)
{
  CH_TIME("MFPoissonOp::residual");
  applyOp(a_lhs, a_phi, a_homogeneous);
  incr(a_lhs, a_rhs, -1);
  scale(a_lhs, -1.0);
}

void MFPoissonOp::computeBoundaryN(const LevelData<MFCellFAB>& a_phi,
                                   bool a_homogeneous)
{
  CH_TIME("MFPoissonOp::computeBoundaryN");
  CH_assert(m_phases == 2);

  LevelData<MFCellFAB>& phi = const_cast<LevelData<MFCellFAB>&>(a_phi);
  phi.exchange(phi.interval());

  DataIterator dit = a_phi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const MFCellFAB& phi = a_phi[dit];
      const EBCellFAB& phi0 = phi.getPhase(0);
      const EBCellFAB& phi1 = phi.getPhase(1);
      const Box& box = phi0.getRegion();
      const EBGraph& ebgraph = phi0.getEBISBox().getEBGraph();
      IntVectSet ivs = ebgraph.getIrregCells(box);
      if (!ivs.isEmpty())
        {
          if (m_ncomp == 1)
            {
              m_jump.computeScalarBV(phi0, phi1,
                                     m_bCoef, dit(),
                                     m_boundaryD[0][dit], m_boundaryD[1][dit],
                                     m_boundaryN[0][dit], m_boundaryN[1][dit],
                                     a_homogeneous);
            }
          else
            {
              CH_assert(m_ncomp==CH_SPACEDIM);
              m_jump.computeVectorBV(phi0, phi1,
                                     m_bCoef, dit(),
                                     m_boundaryD[0][dit], m_boundaryD[1][dit],
                                     m_boundaryN[0][dit], m_boundaryN[1][dit],
                                     a_homogeneous);
            }
        }
    }
}

void MFPoissonOp::setVal(LevelData<MFCellFAB>& a_phi, const Vector<Real> a_values) const
{
  CH_TIME("MFPoissonOp::setVal");
  DataIterator dit = a_phi.dataIterator();
  for (int i=0; i<m_phases; i++)
    {
      aliasMF(*m_alias[0], i, a_phi);
      for (dit.begin(); dit.ok(); ++dit)
        {
          m_alias[0]->operator[](dit).setVal(a_values[i]);
        }
    }
}

void MFPoissonOp::preCond( LevelData<MFCellFAB>&       a_correction,
                           const LevelData<MFCellFAB>& a_residual)
{
  CH_TIME("MFPoissonOp::preCond");
  relax(a_correction, a_residual, 40);
}

void MFPoissonOp::levelJacobi( LevelData<MFCellFAB>&       a_phi,
                               const LevelData<MFCellFAB>& a_rhs)
{
  CH_TIME("MFPoissonOp::levelJacobi");
  LevelData<MFCellFAB>& resid = m_tmp;
  bool homogeneous = true;
  residual(resid, a_phi, a_rhs, homogeneous);

  for (int i=0; i<m_phases; i++)
    {
      aliasMF(*m_alias[0], i, m_weights);
      aliasMF(*m_alias[1], i, resid);
      m_ebops[i]->getInvDiagRHS(*m_alias[0], *m_alias[1]);
    }
  incr(a_phi, m_weights, 0.5);
}

void MFPoissonOp::levelGSRB(LevelData<MFCellFAB>&       a_phi,
                            const LevelData<MFCellFAB>& a_rhs)
{
  CH_TIME("MFPoissonOp::levelGSRB");
  LevelData<MFCellFAB>& resid = m_tmp;

  // Loop over all possibilities
  for (int icolor=0; icolor<=1; ++icolor)
    {
      // Get the residual everywhere
      residual(resid, a_phi, a_rhs, true);
      for (int i=0; i<m_phases; i++)
        {
          aliasMF(*m_alias[0], i, a_phi);
          aliasMF(*m_alias[1], i, resid);
          m_ebops[i]->levelGSRB(*m_alias[0], *m_alias[1], icolor);
        }
    }
}

void MFPoissonOp::levelMulticolorGS( LevelData<MFCellFAB>& a_phi,
                               const LevelData<MFCellFAB>& a_rhs)
{
  CH_TIME("MFPoissonOp::levelMulticolorGS");
  LevelData<MFCellFAB>& resid = m_tmp;

  // Loop over all possibilities (in all dimensions)
  for (int icolor=0; icolor < EBAMRPO_NUMSTEN; ++icolor)
    {
      // Get the residual everywhere
      residual(resid, a_phi, a_rhs, true);
      for (int i=0; i<m_phases; i++)
        {
          aliasMF(*m_alias[0], i, a_phi);
          aliasMF(*m_alias[1], i, resid);
          m_ebops[i]->levelMultiColorGS(*m_alias[0], *m_alias[1], m_colors[icolor]);
        }
    }
}

void MFPoissonOp::applyOp(  LevelData<MFCellFAB>&        a_lhs,
                            const LevelData<MFCellFAB>&  a_phi,
                            bool                         a_homogeneous)
{
  CH_TIME("MFPoissonOp::applyOp");
  computeBoundaryN(a_phi, a_homogeneous);

  for (int i=0; i<m_phases; i++)
  {
    aliasMF(*m_alias[0], i, a_lhs);
    aliasMF(*m_alias[1], i, a_phi);
    m_ebops[i]->applyOp(*m_alias[0], *m_alias[1], NULL,
                        a_homogeneous, true, m_boundaryN+i);
  }
}

void MFPoissonOp::applyOp(  LevelData<MFCellFAB>&        a_lhs,
                            const LevelData<MFCellFAB>&  a_phi,
                            DataIterator&                a_dit,
                            bool                         a_homogeneous)
{
  CH_TIME("MFPoissonOp::applyOp");
  computeBoundaryN(a_phi, a_homogeneous);

  for (int i=0; i<m_phases; i++)
  {
    aliasMF(*m_alias[0], i, a_lhs);
    aliasMF(*m_alias[1], i, a_phi);
    m_ebops[i]->applyOp(*m_alias[0], *m_alias[1], NULL,
                        a_dit, a_homogeneous, true,
                        m_boundaryN+i);
  }
}

void MFPoissonOp::getFlux(MFFluxFAB&                   a_flux,
                          const LevelData<MFCellFAB>&  a_data,
                          const Box&                   a_grid,
                          const DataIndex&             a_dit,
                          Real                         a_scale)
{
  CH_TIME("MFPoissonOp::getFlux");
  for (int iphase=0; iphase<m_phases; iphase++)
    {
      LevelData<EBCellFAB> ebData;
      aliasMF(ebData, iphase, a_data);

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          Box ghostedBox = a_grid;
          ghostedBox.grow(1);
          ghostedBox.grow(idir,-1);
          ghostedBox &= m_ebops[iphase]->getEBLG().getDomain();

          m_ebops[iphase]->getFlux(a_flux.getPhase(iphase),
                                   ebData,
                                   a_grid,
                                   a_dit,
                                   a_scale);

          /*
          // Put the boundary fluxes into a_flux
          const EBGraph& ebgraph = ebData[a_dit].getEBISBox().getEBGraph();
          IntVectSet ivsIrreg = ebgraph.getIrregCells(ebData[a_dit].getRegion());
          VoFIterator vofit(ivsIrreg, ebgraph);
          for (vofit.reset(); vofit.ok(); ++vofit)
            {
              a_flux.getPhase(iphase)[idir]
            }
          */
        }
    }
}

void MFPoissonOp::create(   LevelData<MFCellFAB>& a_lhs,
                            const LevelData<MFCellFAB>& a_rhs)
{
  CH_TIME("MFPoissonOp::create");
  m_ops.create(a_lhs, a_rhs);
}

void MFPoissonOp::createCoarsened(   LevelData<MFCellFAB>& a_lhs,
                                     const LevelData<MFCellFAB>& a_rhs,
                                     const int& a_refRat)
{
  CH_TIME("MFPoissonOp::createCoarsened");
  IntVect ghostVect = a_rhs.ghostVect();
  Vector<EBLevelGrid> eblg(m_phases);
  for (int i=0; i<m_phases; i++)
    {
      eblg[i] = m_ebops[i]->getEBLG();
    }

  DisjointBoxLayout dblCoarsenedFine;
  coarsen(dblCoarsenedFine, eblg[0].getDBL(), a_refRat);
  ProblemDomain coarDom = coarsen(eblg[0].getDomain(), a_refRat);
  IntVect ghostVec = a_rhs.ghostVect();
  Vector<EBISLayout> ebislv(m_phases);
  for (int i=0; i<m_phases; i++)
    {
      eblg[i].getEBIS()->fillEBISLayout(ebislv[i], dblCoarsenedFine, coarDom , ghostVec[0]);
      if (a_refRat > 2)
        {
          ebislv[i].setMaxRefinementRatio(a_refRat, eblg[i].getEBIS());
        }
    }
  //create coarsened data
  Vector<int> vncomp(m_phases, a_rhs.nComp());
  MFCellFactory fact(ebislv, vncomp);
  a_lhs.define(dblCoarsenedFine, a_rhs.nComp(), ghostVec, fact);
}

void MFPoissonOp::assign(   LevelData<MFCellFAB>& a_lhs,
                            const LevelData<MFCellFAB>& a_rhs)
{
  m_ops.assign(a_lhs, a_rhs);
}


Real MFPoissonOp::dotProduct(const LevelData<MFCellFAB>& a_data1,
                             const LevelData<MFCellFAB>& a_data2)
{
  Real accum = 0.0;

  Real volume = 0.0;

  for (DataIterator dit = a_data1.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      Real phaseVolume;
      for (int i=0; i<m_phases; i++)
      {
        const EBCellFAB& data1 = a_data1[d].getPhase(i);
        const EBCellFAB& data2 = a_data2[d].getPhase(i);
        const Box& box = a_data1.getBoxes().get(d);


        accum += EBLevelDataOps::sumKappaDotProduct(phaseVolume,data1,data2,box,
                                    EBLEVELDATAOPS_ALLVOFS,m_domain);

        volume += phaseVolume;
      }
    }



#ifdef CH_MPI
  Real recv;
  int result;

  result = MPI_Allreduce(&accum, &recv, 1, MPI_CH_REAL,
                         MPI_SUM, Chombo_MPI::comm);
  accum = recv;

  result = MPI_Allreduce(&volume, &recv, 1, MPI_CH_REAL,
                         MPI_SUM, Chombo_MPI::comm);
  volume = recv;
#endif

   if (volume > 0.0)
    {
      accum = accum / volume;
    }
  else
    {
      accum = 0.0;
    }

  return accum;
}

void MFPoissonOp::incr( LevelData<MFCellFAB>& a_lhs,
                        const LevelData<MFCellFAB>& a_rhs,
                        Real a_scale)
{
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      a_lhs[d].plus(a_rhs[d], a_scale);
    }
}

void MFPoissonOp::axby( LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_x,
                        const LevelData<MFCellFAB>& a_y,
                        Real a, Real b)
{
  m_ops.axby(a_lhs, a_x, a_y, a, b);
}

void MFPoissonOp::scale(LevelData<MFCellFAB>& a_lhs, const Real& a_scale)
{
  m_ops.scale(a_lhs, a_scale);
}

Real MFPoissonOp::norm(const LevelData<MFCellFAB>& a_x, int a_ord)
{
  Real volume;
  Real rtn = kappaNorm(volume, a_x, a_ord);
  return rtn;
}

Real MFPoissonOp::kappaNorm(Real&                       a_volume,
                            const LevelData<MFCellFAB>& a_data,
                            int                         a_p) const
{

  Real accum = 0.0;

  a_volume = 0.0;
  int ncomp=a_data.nComp();
  for (DataIterator dit = a_data.dataIterator(); dit.ok(); ++dit)
    {
      DataIndex d = dit();
      Vector<Real> cur(ncomp,0.0);
      Real curVolume = 0.0;
      Vector<Real> phaseCur(ncomp,0.0);
      Real phaseVolume;
      for (int i=0; i<m_phases; i++)
      {
        const EBCellFAB& data = a_data[d].getPhase(i);
        const Box& box = a_data.getBoxes().get(d);
        phaseCur = EBLevelDataOps::vectorSumKappaPow(phaseVolume,data,box,
                                                     EBLEVELDATAOPS_ALLVOFS,m_domain,a_p);
        for (int i=0; i<ncomp; i++)
        {
          if (a_p == 0)
          {
            if ( phaseCur[i] > cur[i])
              { cur[i] = phaseCur[i];}
          }

          else
            { cur[i]+= phaseCur[i];
            }
        }
        curVolume += phaseVolume;
      }
      a_volume += curVolume;

      for (int i=0; i<ncomp; i++)
      {
        if (a_p == 0)
          {
            if (cur[i] > accum)
            {
              accum = cur[i];
            }
        }
      else
        {
          accum += cur[i];
        }
      }
    }
#ifdef CH_MPI
  if (a_p == 0)
  {
    Real recv;
    int result;

    result = MPI_Allreduce(&accum, &recv, 1, MPI_CH_REAL,
                           MPI_MAX, Chombo_MPI::comm);
    accum = recv;
  }
  else
  {
    Real recv;
    int result;

    result = MPI_Allreduce(&accum, &recv, 1, MPI_CH_REAL,
                           MPI_SUM, Chombo_MPI::comm);
    accum = recv;

    result = MPI_Allreduce(&a_volume, &recv, 1, MPI_CH_REAL,
                           MPI_SUM, Chombo_MPI::comm);
    a_volume = recv;
  }
#endif
  if (a_p != 0)
    {
      if (a_volume > 0.0)
        {
          accum = accum / a_volume;
        }
      else
        {
          accum = 0;
        }
      Real invPow = 1.0/a_p;
      accum = pow(accum,invPow);
    }

  return accum;
}


void MFPoissonOp::setToZero( LevelData<MFCellFAB>& a_x)
{
  m_ops.setToZero(a_x);
}

/*@}*/

/**
   \name MGLevelOp functions */
/*@{*/

void MFPoissonOp::relax(LevelData<MFCellFAB>& a_e,
                        const LevelData<MFCellFAB>& a_residual,
                        int iterations)
{
  for (int i=0; i<iterations; i++)
    {
      if (m_relax == 0)
        {
          levelJacobi(a_e, a_residual);
        }
      else if (m_relax == 1)
        {
          levelMulticolorGS(a_e, a_residual);
        }
      else if (m_relax == 2)
        {
          levelGSRB(a_e, a_residual);
        }
      else
        {
          MayDay::Error(" MFPoissoOp: Invalid relaxation type");
        }
    }
}

void MFPoissonOp::createCoarser(LevelData<MFCellFAB>& a_coarse,
                                const LevelData<MFCellFAB>& a_fine,
                                bool ghosted)
{
  CH_TIME("MFPoissonOp::createCoarser");
  Vector<EBISLayout> ebislv(m_phases);
  for (int i=0; i<m_phases; i++)
    {
      ebislv[i] = m_ebops[i]->getEBLGCoarMG().getEBISL();
    }
  Vector<int> vncomp(m_phases, m_ncomp);
  MFCellFactory factory(ebislv, vncomp);
  DisjointBoxLayout dbllv = m_ebops[0]->getEBLGCoarMG().getDBL();
  a_coarse.define(dbllv, m_ncomp, a_fine.ghostVect(), factory);
}

/**
   calculate restricted residual
   a_resCoarse[2h] = I[h->2h] (rhsFine[h] - L[h](phiFine[h])
*/
void MFPoissonOp::restrictResidual(LevelData<MFCellFAB>& a_resCoarse,
                                   LevelData<MFCellFAB>& a_phiFine,
                                   const LevelData<MFCellFAB>& a_rhsFine)
{
  CH_TIME("MFPoissonOp::restrictResidual");
  for (int i=0; i<m_phases; i++)
  {
    aliasMF(*m_alias[0], i, a_resCoarse);
    aliasMF(*m_alias[1], i, a_phiFine);
    aliasMF(*m_alias[2], i, a_rhsFine);
    m_ebops[i]->restrictResidual(*m_alias[0], *m_alias[1], *m_alias[2]);
  }
}

/**
   correct the fine solution based on coarse correction
   a_phiThisLevel += I[2h->h](a_correctCoarse)
*/
void MFPoissonOp::prolongIncrement(LevelData<MFCellFAB>& a_phiThisLevel,
                                   const LevelData<MFCellFAB>& a_correctCoarse)
{
  CH_TIME("MFPoissonOp::prolongIncrement");
  for (int i=0; i<m_phases; i++)
  {
    aliasMF(*m_alias[0], i, a_phiThisLevel);
    aliasMF(*m_alias[1], i, a_correctCoarse);
    m_ebops[i]->prolongIncrement(*m_alias[0], *m_alias[1]);
  }
}
/*@}*/

/**
   \name AMRLevelOp functions */
/*@{*/

/** a_residual = a_rhs - L(a_phi, a_phiFine, a_phiCoarse) */
void MFPoissonOp::AMRResidual(LevelData<MFCellFAB>& a_residual,
                              const LevelData<MFCellFAB>& a_phiFine,
                              const LevelData<MFCellFAB>& a_phi,
                              const LevelData<MFCellFAB>& a_phiCoarse,
                              const LevelData<MFCellFAB>& a_rhs,
                              bool a_homogeneousBC,
                              AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp)
{
  CH_TIME("MFPoissonOp::AMRResidual");

  AMROperator(a_residual, a_phiFine, a_phi, a_phiCoarse,
              a_homogeneousBC, a_finerOp);

  scale(a_residual,-1.0);
  incr(a_residual,a_rhs,1.0);

}
/** residual assuming no more coarser AMR levels */

void MFPoissonOp::AMRResidualNC(LevelData<MFCellFAB>& a_residual,
                                const LevelData<MFCellFAB>& a_phiFine,
                                const LevelData<MFCellFAB>& a_phi,
                                const LevelData<MFCellFAB>& a_rhs,
                                bool a_homogeneousBC,
                                AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp)
{
  CH_TIME("MFPoissonOp::AMRResidualNC");

  AMROperatorNC(a_residual, a_phiFine, a_phi, a_homogeneousBC,
                a_finerOp);

  scale(a_residual,-1.0);
  incr(a_residual,a_rhs,1.0);
}

/** a_residual = a_rhs - L(a_phi, a_phiCoarse)  */
void MFPoissonOp::AMRResidualNF(LevelData<MFCellFAB>& a_residual,
                                const LevelData<MFCellFAB>& a_phi,
                                const LevelData<MFCellFAB>& a_phiCoarse,
                                const LevelData<MFCellFAB>& a_rhs, bool a_homogeneousBC)
{
  CH_TIME("MFPoissonOp::AMRResidualNF");

  AMROperatorNF(a_residual, a_phi, a_phiCoarse, a_homogeneousBC);

  scale(a_residual,-1.0);
  incr(a_residual,a_rhs,1.0);
}



/** a_operator = a_rhs - L(a_phi, a_phiFine, a_phiCoarse) */
void MFPoissonOp::AMROperator(LevelData<MFCellFAB>& a_LofPhi,
                              const LevelData<MFCellFAB>& a_phiFine,
                              const LevelData<MFCellFAB>& a_phi,
                              const LevelData<MFCellFAB>& a_phiCoarse,
                              bool a_homogeneousBC,
                              AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp)
{
  CH_TIME("MFPoissonOp::AMROperator");
  computeBoundaryN(a_phi, a_homogeneousBC);
  for (int i=0; i<m_phases; i++)
    {
      aliasMF(*m_alias[0], i, a_LofPhi);
      aliasMF(*m_alias[1], i, a_phi);
      aliasMF(*m_alias[2], i, a_phiCoarse);
      aliasMF(*m_alias[3], i, a_phiFine);
      m_ebops[i]->applyOp(*m_alias[0], *m_alias[1], m_alias[2],
                          a_homogeneousBC, false,
                          m_boundaryN+i);
      MFPoissonOp* finerOp = (MFPoissonOp*)a_finerOp;
      m_ebops[i]->reflux(*m_alias[0], *m_alias[3], *m_alias[1],
                         finerOp->m_ebops[i]);
    }
}

/** AMR operator assuming no more coarser AMR levels */

void MFPoissonOp::AMROperatorNC(LevelData<MFCellFAB>& a_LofPhi,
                                const LevelData<MFCellFAB>& a_phiFine,
                                const LevelData<MFCellFAB>& a_phi,
                                bool a_homogeneousBC,
                                AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp)
{
  CH_TIME("MFPoissonOp::AMROperatorNC");
  computeBoundaryN(a_phi, a_homogeneousBC);
  for (int i=0; i<m_phases; i++)
    {
      aliasMF(*m_alias[0], i, a_LofPhi);
      aliasMF(*m_alias[1], i, a_phi);
      aliasMF(*m_alias[3], i, a_phiFine);
      m_ebops[i]->applyOp(*m_alias[0], *m_alias[1], NULL,
                          a_homogeneousBC, true,
                          m_boundaryN+i);
      MFPoissonOp* finerOp = (MFPoissonOp*)a_finerOp;
      m_ebops[i]->reflux(*m_alias[0], *m_alias[3], *m_alias[1],
                         finerOp->m_ebops[i]);
    }
}

/** a_residual = a_rhs - L(a_phi, a_phiCoarse)  */
void MFPoissonOp::AMROperatorNF(LevelData<MFCellFAB>& a_LofPhi,
                                const LevelData<MFCellFAB>& a_phi,
                                const LevelData<MFCellFAB>& a_phiCoarse,
                                bool a_homogeneousBC)
{
  CH_TIME("MFPoissonOp::AMROperatorNF");
  computeBoundaryN(a_phi, a_homogeneousBC);
  for (int i=0; i<m_phases; i++)
    {
      aliasMF(*m_alias[0], i, a_LofPhi);
      aliasMF(*m_alias[1], i, a_phi);
      aliasMF(*m_alias[2], i, a_phiCoarse);
      m_ebops[i]->applyOp(*m_alias[0], *m_alias[1], m_alias[2],
                          a_homogeneousBC, false,
                          m_boundaryN+i);
    }
}

/** a_resCoarse = I[h-2h]( a_residual - L(a_correction, a_coarseCorrection))
    it is assumed that a_resCoarse has already been filled in with the coarse
    version of AMRResidualNF and that this operation is free to overwrite
    in the overlap regions.
*/

void MFPoissonOp::AMRRestrict(LevelData<MFCellFAB>& a_resCoarse,
                              const LevelData<MFCellFAB>& a_residual,
                              const LevelData<MFCellFAB>& a_correction,
                              const LevelData<MFCellFAB>& a_coarseCorrection,
                              bool                        a_skip_res)
{
  CH_TIME("MFPoissonOp::AMRRestrict");
  for (int i=0; i<m_phases; i++)
  {
    aliasMF(*m_alias[0], i, a_resCoarse);
    aliasMF(*m_alias[1], i, a_residual);
    aliasMF(*m_alias[2], i, a_correction);
    aliasMF(*m_alias[3], i, a_coarseCorrection);
    m_ebops[i]->AMRRestrict(*m_alias[0], *m_alias[1], *m_alias[2], *m_alias[3]);
  }

}

/** a_correction += I[h->h](a_coarseCorrection) */
void MFPoissonOp::AMRProlong(LevelData<MFCellFAB>& a_correction,
                             const LevelData<MFCellFAB>& a_coarseCorrection)
{
  CH_TIME("MFPoissonOp::AMRProlong");
  for (int i=0; i<m_phases; i++)
  {
    aliasMF(*m_alias[0], i, a_correction);
    aliasMF(*m_alias[1], i, a_coarseCorrection);
    m_ebops[i]->AMRProlong(*m_alias[0], *m_alias[1]);
  }
}

/** a_residual = a_residual - L(a_correction, a_coarseCorrection) */
void MFPoissonOp::AMRUpdateResidual(LevelData<MFCellFAB>& a_residual,
                                    const LevelData<MFCellFAB>& a_correction,
                                    const LevelData<MFCellFAB>& a_coarseCorrection)
{
  CH_TIME("MFPoissonOp::AMRUpdateResidual");
  bool homogeneousBC = true;
  computeBoundaryN(a_correction, homogeneousBC);
  for (int i=0; i<m_phases; i++)
  {
    aliasMF(*m_alias[0], i, a_residual);
    aliasMF(*m_alias[1], i, a_correction);
    aliasMF(*m_alias[2], i, a_coarseCorrection);
    m_ebops[i]->AMRUpdateResidual(*m_alias[0], *m_alias[1], *m_alias[2], m_boundaryN+i);
  }
}

///
/**
   compute norm over all cells on coarse not covered by finer
*/
Real MFPoissonOp::AMRNorm(const LevelData<MFCellFAB>& a_coarseResid,
                          const LevelData<MFCellFAB>& a_fineResid,
                          const int&                   a_refRat,
                          const int&                   a_ord)
{
  CH_TIME("MFPoissonOp::AMRNorm");
  Real m = 0;
  for (int i=0; i<m_phases; i++)
  {
    aliasMF(*m_alias[0], i, a_coarseResid);
    aliasMF(*m_alias[1], i, a_fineResid);
    Real norm = m_ebops[i]->AMRNorm(*m_alias[0], *m_alias[1], a_refRat, a_ord);

    m = Max(m, norm);

  }
  return m;
}

#include "NamespaceFooter.H"

#ifndef ____Debugging_Tools____
#define ____Debugging_Tools____
#define DBG_LOX  -2
#define DBG_HIX  3
#define DBG_LOY  12
#define DBG_HIY  17
#define DBG_LOZ  0
#define DBG_HIZ  0
#define VALID_ONLY false

/**
//   for dumpFracPr and dumpFracCh
//   Makes box = defined fraction of the domain.
// set the numerators and denominators of fraction you want displayed
// This:
//#define DBG_NUMER_LO 1
//#define DBG_DENOM_LO 8
//#define DBG_NUMER_HI 1
//#define DBG_DENOM_HI 4
//  Will display the fraction of the domain that is between 1/8 and 1/4 of
//   of the domain area.  (so near the lower left corner)
//  This may go down in history as the worst API ever but it is helpful when you
//  need to keep an eye on a particular region of the domain when AMR is involved.
**/
#define DBG_NUMER_LO_X 0
#define DBG_DENOM_LO_X 1
#define DBG_NUMER_HI_X 1
#define DBG_DENOM_HI_X 8

#define DBG_NUMER_LO_Y 3
#define DBG_DENOM_LO_Y 16
#define DBG_NUMER_HI_Y 1
#define DBG_DENOM_HI_Y 2

#define DBG_NUMER_LO_Z 0
#define DBG_DENOM_LO_Z 1
#define DBG_NUMER_HI_Z 0
#define DBG_DENOM_HI_Z 1


#include "FArrayBox.H"
#include "LevelData.H"
#include "DisjointBoxLayout.H"



/// Old school Debugging utilities 
/**
 // If you want to print areas of data to screen during a debug session,
 // his is the way you do it.  All of these functions can be called from gdb.
 // TLDR:
 // This is a set  of crude but effective tools to  print out stuff from gdb.   
 // Send pointers to the data to these functions and the data gets printed out neatly.    
 // The user is meant to hack the the bounds according to her needs.   
 // The macros above are the primary buttons to modify.    
 // 
 // Set the bounds (DBG_LOX and so on) to customize which area gets printed.
 // Set VALID_ONLY to TRUE if you only want valid data printed when a 
 //     level gets printed out (sometimes it is help ful to see ghost data, other times not so much).
 // 
 // This command:
 // (gdb) p PrhDebug::dumpLDAreaPr(&a_cor)  
 // yeilds:
 // ibox = 0:
 // data_box = [(-2,-2),( 3, 3)], area_box = [( 0, 0),( 3, 3)]
 // data for variable 0:
 // ( 0, 3): 4.86993301e-05  ( 1, 3): 3.13317884e-04  ( 2, 3): 1.72362925e-04  ( 3, 3): 6.44824962e-05  
 // ( 0, 2): 3.13287135e-04  ( 1, 2): 5.52541534e-04  ( 2, 2): 8.82170603e-04  ( 3, 2): 1.72362925e-04  
 // ( 0, 1): 1.72411307e-04  ( 1, 1): 8.82090868e-04  ( 2, 1): 5.52541534e-04  ( 3, 1): 3.13317884e-04  
 // ( 0, 0): 6.46298245e-05  ( 1, 0): 1.72411307e-04  ( 2, 0): 3.13287135e-04  ( 3, 0): 4.86993301e-05  
 // and so on.
 **/
namespace PrChDebug
{
  typedef   Chombo::IntVect                      ch_iv;
  typedef   Chombo::Box                          ch_box;
  typedef   Chombo::FArrayBox                    ch_fab;
  typedef   Chombo::DataIterator                 ch_dit;
  typedef   Chombo::DisjointBoxLayout            ch_dbl;
  typedef   Chombo::LevelData<ch_fab>            ch_ldf;
                                                  
  
  void
  dumpSten(pr_sten* a_sten_ptr)
  {
    if(a_sten_ptr != NULL)
    {
      a_sten_ptr->print();
    }
  }

  void 
  dumpCh(ch_fab* a_dataPtr, ch_box a_areaBox)
  {

    if(a_dataPtr != NULL)
    {
      using std::cout;
      using std::end;
      cout    << setprecision(8)
              << setiosflags(ios::showpoint)
              << setiosflags(ios::scientific);

      ch_box area = a_areaBox;

      ch_fab & data    =*a_dataPtr;
      ch_box databox   = a_dataPtr->box();
      ch_iv  hidata    =   databox.bigEnd();
      ch_iv  lodata    =   databox.smallEnd();

      ///to get around weird box types
      ch_box printBox(lodata, hidata);
      printBox &= a_areaBox;
      ch_iv lodeb = printBox.smallEnd();
      ch_iv hideb = printBox.bigEnd();
      
      for(int ivar = 0; ivar < data.nComp(); ivar++)
      {
        cout << "data_box = " << databox << ", area_box = " << a_areaBox <<  ", ivar = " << ivar << endl;
#if DIM==3        
        for(int k = hideb[2]; k >= lodeb[2]; k--)
        {
#endif          
          for(int j = hideb[1]; j >= lodeb[1]; j--)
          {
            for(int i = lodeb[0]; i <= hideb[0]; i++)
            {
              ch_iv pt(D_DECL(i,j,k));
              cout << pt << ":";
              double data_value = data(pt, ivar);
              //keeps positive and negative numbers lined up
              if(    data_value >= 0)
              {
                cout << " "; 
              }
              // 
              cout  << data_value << "  ";
            } //end loop in x dir
            cout << endl;
          } //end loop in  y dir
#if DIM==3
        } //end loop in z direction
#endif
      }// end loop over components
    } //end if not null
  } //end function dumpAreaCh

  ///
  void 
  dumpAreaCh(ch_fab* dataPtr)
  {
    static const ch_iv       s_areaLo(D_DECL( DBG_LOX, DBG_LOY, DBG_LOZ));
    static const ch_iv       s_areaHi(D_DECL( DBG_HIX, DBG_HIY, DBG_HIZ));
    static const ch_box     s_areaBox(s_areaLo, s_areaHi);
    dumpCh(dataPtr, s_areaBox);
  }
  
  
  /// print to screen (very neatly) data from the area whereever it is found in chombo data
  void 
  dumpLDAreaCh(ch_ldf* dataPtr)
  {
    if(dataPtr != NULL)
    {
      static const ch_iv       s_areaLo(D_DECL( DBG_LOX, DBG_LOY, DBG_LOZ));
      static const ch_iv       s_areaHi(D_DECL( DBG_HIX, DBG_HIY, DBG_HIZ));
      static const ch_box     s_areaBox(s_areaLo, s_areaHi);

      ch_ldf & data = *dataPtr;
      ch_dbl   grid = data.disjointBoxLayout();
      ch_dit  dit = data.dataIterator();
      for(int ibox = 0; ibox < dit.size(); ibox++)
      {
        ch_box valid = grid[dit[ibox]];
        ch_box valid_area = valid & s_areaBox;
        cout << "ibox = " << ibox << ", valid = " << valid << ", valid_area = :" << valid_area << endl;
        bool valid_only = VALID_ONLY;

        if(valid_only)
        {
          dumpCh(&data[dit[ibox]], valid_area);
        }
        else if(!valid_area.isEmpty())
        {
          dumpCh(&data[dit[ibox]], s_areaBox);
        }
      }
    }
  }
  void dumpLayoutCh(ch_dbl* a_data_ptr)
  {
    if(a_data_ptr != NULL)
    {
      a_data_ptr->print();
    }
  }
  void dumpLDLayoutCh(ch_ldf* a_data_ptr)
  {
    if(a_data_ptr != NULL)
    {
      a_data_ptr->disjointBoxLayout().print();
    }
  }
}//end namespace PrChDebug

///
inline void dummy_function_to_get_functions_into_symbol_table()
{
  PrChDebug::dumpLayoutCh(NULL);
  PrChDebug::dumpLDLayoutCh(NULL);
  PrChDebug::dumpAreaCh(NULL);
  PrChDebug::dumpLDAreaCh(NULL);
  PrChDebug::dumpSten(NULL);
}

#undef DBG_LOX 
#undef DBG_LOY 
#undef DBG_LOZ 
#undef DBG_HIX 
#undef DBG_HIY 
#undef DBG_HIZ 
#undef VALID_ONLY
#undef DBG_NUMER_LO 
#undef DBG_DENOM_LO 
#undef DBG_NUMER_HI 
#undef DBG_DENOM_HI 


#endif

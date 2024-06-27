/*--------------------------------------------------------------------------*/
/*------------------------- File MCFLemonSolver.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the MCFLemonSolver class.
 *
 * \author Daniele Caliandro \n
 *         Universita' di Pisa \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Daniele Caliandro, Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

<<<<<<< HEAD
=======
//#include <math.h>

>>>>>>> 3587a0f95544d63e767d768c535e345b48386907
#include "MCFLemonSolver.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/
using namespace SMSpp_di_unipi_it;


<<<<<<< HEAD
// register the various LEMONSolver GR , V , C > to the Solver factory

SMSpp_insert_in_factory_cpp_0_t(
    MCFLemonSolverNetworkSimplex<SmartDigraph, double, double>);

SMSpp_insert_in_factory_cpp_0_t(
    MCFLemonSolverNetworkSimplex<SmartDigraph, long, long>);

SMSpp_insert_in_factory_cpp_0_t(
    MCFLemonSolverNetworkSimplex<MCFListDigraph, double, double>);

SMSpp_insert_in_factory_cpp_0_t(
    MCFLemonSolverNetworkSimplex<MCFListDigraph, long, long>);

SMSpp_insert_in_factory_cpp_0_t(
    MCFLemonSolverCycleCanceling<SmartDigraph, double, double>);

SMSpp_insert_in_factory_cpp_0_t(
    MCFLemonSolverCycleCanceling<SmartDigraph, long, long>);
=======
template<typename GR, typename V, typename C>
struct Fields< SMSppCapacityScaling< GR, V, C> > {
  int f_factor;
};

template<typename GR, typename V, typename C>  
struct Fields< SMSppCostScaling< GR, V, C> > {
  using CSMethod = typename SMSppCostScaling< GR, V, C >::Method;
  int f_factor;
  CSMethod f_method;
};
template<typename GR, typename V, typename C>
struct Fields< NetworkSimplex<GR, V, C> > {
  using NSPivotRule = typename NetworkSimplex<GR, V, C>::PivotRule;
  NSPivotRule *f_pivot_rule;
};
template<typename GR, typename V, typename C>
struct Fields< CycleCanceling<GR, V, C> > {
  using CCMethod = typename CycleCanceling<GR, V, C>::Method;
  CCMethod *f_method;
};



>>>>>>> 3587a0f95544d63e767d768c535e345b48386907

SMSpp_insert_in_factory_cpp_0_t(
    MCFLemonSolverCycleCanceling<MCFListDigraph, double, double>);

SMSpp_insert_in_factory_cpp_0_t(
    MCFLemonSolverCycleCanceling<MCFListDigraph, long, long>);

<<<<<<< HEAD
SMSpp_insert_in_factory_cpp_0_t(
    MCFLemonSolverCostScaling<SmartDigraph, long, long>);
=======
// SMSpp_insert_in_factory_cpp_0( MCFSolverState );
>>>>>>> 3587a0f95544d63e767d768c535e345b48386907

SMSpp_insert_in_factory_cpp_0_t(
    MCFLemonSolverCostScaling<SmartDigraph, double, double>);

SMSpp_insert_in_factory_cpp_0_t(
    MCFLemonSolverCostScaling<MCFListDigraph, double, double>);

SMSpp_insert_in_factory_cpp_0_t(
    MCFLemonSolverCostScaling<MCFListDigraph, long, long>);

SMSpp_insert_in_factory_cpp_0_t(
    MCFLemonSolverCapacityScaling<SmartDigraph, double, double>);

SMSpp_insert_in_factory_cpp_0_t(
    MCFLemonSolverCapacityScaling<SmartDigraph, long, long>);

SMSpp_insert_in_factory_cpp_0_t(
    MCFLemonSolverCapacityScaling<SmartDigraph, long, double>);

SMSpp_insert_in_factory_cpp_0_t(
    MCFLemonSolverCapacityScaling<MCFListDigraph, double, double>);

SMSpp_insert_in_factory_cpp_0_t(
    MCFLemonSolverCapacityScaling<MCFListDigraph, long, long>);

SMSpp_insert_in_factory_cpp_0_t(
    MCFLemonSolverCapacityScaling<MCFListDigraph, long, double>);



/*--------------------------------------------------------------------------*/
/*-------------------------------- SPECIALIZED CLASSES --------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------------- NETWORKSIMPLEX --------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/


/** Specialized MCFLemonSolver<NetworkSimplex, GR, V, C> that contains specialized
 *  compute() method, enums for indexing algorithimc parameters and function
 *  set/get_*_par for manage them.
 *  This class derives from MCFLemonSolverBase for using set_Block() methods and
 *  for inherit f_algo and dgp fields.
 * 
 *  Template parameters are:
 * 
 *  - NetworkSimplex, implements the primal Network Simplex algorithm for finding
 *    a minimum cost flow. This algorithm is a highly efficient specialized version
 *    of the linear programming simplex method directly for the minimum cost flow
 *    problem.
 * 
 *  - GR that represents the directed graph, the possibilities are described in the
 *    file MCFLemonSolver.h at line 339
 *  
 *  - V, which is the type of flows / deficits; typically, double can be used
 *    for maximum compatibility, but int (or even smaller) would yeld better
 *    performances;
 *  
 * - C, which is the type of ar costs; typically, double can be used for
 *    maximum compatibility, but int (or even smaller) would yeld better
 *    performances;
*/

template< typename GR, typename V, typename C>
class MCFLemonSolver<NetworkSimplex, GR, V, C> : virtual public CDASolver, public MCFLemonSolverBase<NetworkSimplex, GR, V, C>, Fields< NetworkSimplex<GR, V, C> >
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

   

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
   kUnEval = 0     compute() has not been called yet

   kUnbounded = kUnEval + 1     the model is provably unbounded
 *  @{ */
  const int kErrorStatus = -1;
/*
    
    kUnEval = 0     compute() has not been called yet

    kUnbounded = kUnEval + 1     the model is provably unbounded

    kInfeasible                  the model is provably infeasible

    kBothInfeasible = kInfeasible + 1     both primal and dual infeasible

    kOK = 7         successful compute()
                    Any return value between kUnEval (excluded) and kOK
        (included) means that the object ran smoothly

    kStopTime = kOK + 1          stopped because of time limit

    kStopIter                    stopped because of iteration limit

    kError = 15     compute() stopped because of unrecoverable error
                    Any return value >= kError means that the object was
         forced to stop due to some error, e.g. of numerical nature

    kLowPrecision = kError + 1   a solution found but not provably optimal
    */



/** @} ---------------------------------------------------------------------*/
/*-------------- CONSTRUCTING AND DESTRUCTING MCFLemonSolver ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing MCFLemonSolver
 *  @{ */

 /// constructor: Initializes algorithm parameters
 /** Void constructor. Define f_pivot_rule_type to the default algorithm parameters used by NetworkSimplex.
  */

 MCFLemonSolver( void ) :  CDASolver(), MCFLemonSolverBase<NetworkSimplex, GR, V , C>(), Fields<NetworkSimplex<GR, V, C > >() {
  f_pivot_rule_type = NetworkSimplex<GR, V, C>::PivotRule::BLOCK_SEARCH;
  }
 
  /// destructor: delete algorithm parameter
  /** Does nothing special, delete Fields f_pivot_rule, an algorithmic parameter of
   *  NetworkSimplex
   * */  
  ~MCFLemonSolver( void ) {
    delete Fields<NetworkSimplex<GR, V, C>>::f_pivot_rule;  
  }
 
 

/*--------------------------------------------------------------------------*/
 /* intMaxIter = 0     maximum iterations for the next call to solve()

    intMaxSol          maximum number of different solutions to report

    intLogVerb         "verbosity" of the log

    intMaxDSol         maximum number of different dual solutions

    intLastParCDAS     first allowed parameter value for derived classes
    */

/*--------------------------------------------------------------------------*/
 /* dblMaxTime = 0    maximum time for the next call to solve()

    dblRelAcc         relative accuracy for declaring a solution optimal

    dblAbsAcc          absolute accuracy for declaring a solution optimal

    dblUpCutOff        upper cutoff for stopping the algorithm

    dblLwCutOff        lower cutoff for stopping the algorithm

    dblRAccSol          maximum relative error in any reported solution

    dblAAccSol          maximum absolute error in any reported solution

    dblFAccSol          maximum constraint violation in any reported solution

    dblRAccDSol         maximum relative error in any dual solution

    dblAAccDSol         maximum absolute error in any dual solution

    dblFAccDSol         maximum absolute error in any dual solution

    dblLastParCDAS      first allowed parameter value for derived classes
    */

/*--------------------------------------------------------------------------*/
 enum LEMON_NS_dbl_par_type{
        dblLastParLEMON_NS ///< first allowed parameter value for derived classes
                           /**< convenience value for easily allow derived classes
                            * to further extend the set of types of return codes */
 };

 enum LEMON_NS_int_par_type{
        kPivot = intLastParCDAS,
        intLastParLEMON_NS ///< first allowed parameter value for derived classes
                           /**< convenience value for easily allow derived classes
                            * to further extend the set of types of return codes */
 };
 
 enum str_par_type_LEMON_NS {
  strDMXFile = strLastParCDAS ,  ///< DMX filename to output the instance
  strLastParLEMON_NS    ///< first allowed parameter value for derived classes
                   /**< convenience value for easily allow derived classes
                    * to further extend the set of types of return codes */
  };  



/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *
 * Parameter-wise, MCFSolver maps the parameters of [CDA]Solver
 *
 *  intMaxIter = 0    maximum iterations for the next call to solve()
 *  intMaxSol         maximum number of different solutions to report
 *  intLogVerb        "verbosity" of the log
 *  intMaxDSol        maximum number of different dual solutions
 *
 *  dblMaxTime = 0    maximum time for the next call to solve()
 *  dblRelAcc         relative accuracy for declaring a solution optimal
 *  dblAbsAcc         absolute accuracy for declaring a solution optimal
 *  dblUpCutOff       upper cutoff for stopping the algorithm
 *  dblLwCutOff       lower cutoff for stopping the algorithm
 *  dblRAccSol        maximum relative error in any reported solution
 *  dblAAccSol        maximum absolute error in any reported solution
 *  dblFAccSol        maximum constraint violation in any reported solution
 *  dblRAccDSol       maximum relative error in any dual solution
 *  dblAAccDSol       maximum absolute error in any dual solution
 *  dblFAccDSol       maximum absolute error in any dual solution
 *
 * into the parameter of MCFClass
 *
 * kMaxTime = 0       max time
 * kMaxIter           max number of iteration
 * kEpsFlw            tolerance for flows
 * kEpsDfct           tolerance for deficits
 * kEpsCst            tolerance for costs
 *
 * It then "extends" them, using
 *
 *  intLastParCDAS    first allowed parameter value for derived classes
 *  dblLastParCDAS    first allowed parameter value for derived classes
 *
 * In particular, one now has
 *
 * intLastParCDAS ==> kReopt             whether or not to reoptimize
 *
 * and any other parameter of specific :MCFClass following. This is done
 * via the two const static arrays Solver_2_MCFClass_int and
 * Solver_2_MCFClass_dbl, with a negative entry meaning "there is no such
 * parameter in MCFSolver".
 *
 *  @{ */


    /*--------------------------------------------------------------------------*/
    // set the ostream for the Solver log
    // not really, MCFClass objects are remarkably silent
    //
    // virtual void set_log( std::ostream *log_stream = nullptr ) override;

    /*--------------------------------------------------------------------------*/
 
  /// @brief set the parameter par with value
  /// @param par 
  /// @param value 
  void set_par(idx_type par, int value) override {
        
    if(par == kPivot){      

      if( (value < 0) || (value > 4)){
        throw(std::invalid_argument(std::to_string(value)));
      }                

      if( value == f_pivot_rule_type ){
        return; //nothing is changed
      }
                
      f_pivot_rule_type = value;
      delete Fields<NetworkSimplex<GR, V, C>>::f_pivot_rule;
      Fields<NetworkSimplex<GR, V, C>>::f_pivot_rule = NULL;
      return;
    }

    CDASolver::set_par(par, value);
    return;

  }
    
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    //TODO: change MCFC function to Algo function.
    /*void set_par(idx_type par, double value) override
    {
      if (Solver_2_MCFClass_dbl[par] >= 0)
      ;
      //  MCFC::SetPar(Solver_2_MCFClass_dbl[par], double(value));
    }*/
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/


    /** @} ---------------------------------------------------------------------*/
    /*--------------------- METHODS FOR SOLVING THE Block ----------------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Solving the MCF encoded by the current MCFBlock
     *  @{ */
    /// (try to) solve the MCF encoded in the MCFBlock 
        int compute( bool changedvars = true ) override {
        const static std::array<int, 6> LemonStatus_2_MCFstatus = {
                kErrorStatus, SMSpp_di_unipi_it::LEMON_sol_type::OPTIMAL, kErrorStatus , 
                SMSpp_di_unipi_it::LEMON_sol_type::INFEASIBLE,
                SMSpp_di_unipi_it::LEMON_sol_type::UNBOUNDED, kErrorStatus};
        
        const static std::array<int, 6> MCFstatus_2_sol_type = {
                kUnEval, Solver::kOK, kStopTime, kInfeasible, Solver::kUnbounded,
                Solver::kError};

        lock(); // first of all, acquire self-lock

        if (!f_Block)            // there is no [MCFBlock] to solve
                return (kBlockLocked); // return error

        bool owned = f_Block->is_owned_by(f_id); // check if already locked
        if ((!owned) && (!f_Block->read_lock())) // if not try to read_lock
                return (kBlockLocked);                 // return error on failure

        // while [read_]locked, process any outstanding Modification
        //TODO: ensure that modification are actually processed for MCFLemonSolver.
        //process_outstanding_Modification();
        
        if (!f_dmx_file.empty())
        { // if so required
                // output the current instance (after the changes) to a DMX file
                std::ofstream ProbFile(f_dmx_file, ios_base::out | ios_base::trunc);
                if (!ProbFile.is_open())
                throw(std::logic_error("cannot open DMX file " + f_dmx_file));

                //WriteMCF(ProbFile);
                writeDimacsMat(ProbFile, *MCFLemonSolver::dgp);
                ProbFile.close();
        }

        if (!owned)               // if the [MCF]Block was actually read_locked
                f_Block->read_unlock(); // read_unlock it

        auto start = chrono::system_clock::now();
        
        if(Fields<NetworkSimplex< GR, V, C > >::f_pivot_rule != NULL){
        this->status = MCFLemonSolver::f_algo->run(*Fields<NetworkSimplex< GR, V, C > >::f_pivot_rule);
        }else{
                //Build f_pivot and execute run() method
                Fields<NetworkSimplex<GR, V, C> >::f_pivot_rule = new typename NetworkSimplex<GR, V, C>::PivotRule(static_cast<typename NetworkSimplex<GR, V, C>::PivotRule>(f_pivot_rule_type));
                this->status = MCFLemonSolver::f_algo->run(*Fields<NetworkSimplex< GR, V, C > >::f_pivot_rule);
        }
        
        auto end = chrono::system_clock::now();

        chrono::duration< double > elapsed = end - start;
        ticks = elapsed.count();

        
        
        unlock(); // release self-lock

        // now give out the result: note that the vector MCFstatus_2_sol_type[]
        // starts from 0 whereas the first value of MCFStatus is -1 (= kUnSolved),
        // hence the returned status has to be shifted by + 1
        return (MCFstatus_2_sol_type[this->get_status()]);
        }

    

    /** @} ---------------------------------------------------------------------*/
    /*---------------------- METHODS FOR READING RESULTS -----------------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Accessing the found solutions (if any)
     *  @{ */
    
    int get_status(void) const 
    {
      return (this->status);
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    ///Get elapsed time for run() method
    double get_elapsed_time(void) const override
    {
      return (this->ticks);
    }



    /*--------------------------------------------------------------------------*/
    //Return the lower bound solution(optimal) for the problem
    OFValue get_lb(void) override { return OFValue(MCFLemonSolver::f_algo->totalCost()); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    //Return the upper bound solution(optimal) for the problem
    OFValue get_ub(void) override { return OFValue(MCFLemonSolver::f_algo->totalCost()); }

    /*--------------------------------------------------------------------------*/
    bool has_var_solution(void) override
    {
     switch (this->get_status()) {
      case( NetworkSimplex< GR , V , C >::ProblemType::OPTIMAL ):
      case( NetworkSimplex< GR , V , C >::ProblemType::UNBOUNDED ):  
       return( true );
      default:
       return( false );
      }
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    bool has_dual_solution(void) override
    {
     switch( this->get_status() ) {
      case( NetworkSimplex< GR , V , C >::ProblemType::OPTIMAL ):
      case( NetworkSimplex< GR , V , C >::ProblemType::INFEASIBLE ):
       return( true );
      default:
       return( false );
      }
    }

    /*--------------------------------------------------------------------------*/
    /*
     virtual bool is_var_feasible( void ) override { return( true ); }

     virtual bool is_dual_feasible( void ) override { return( true ); }
    */
    /*--------------------------------------------------------------------------*/
    /// write the "current" flow in the x ColVariable of the MCFBlock
    /** Write the "current" flow in the x ColVariable of the MCFBlock. To keep
     * the same format as MCFBlock::get_Solution() and
     * MCFBlock::map[forward/back]_Modification(), the Configuration *solc can
     * be used to "partly" save it. In particular, if solc != nullptr, it is
     * a SimpleConfiguration< int >, and solc->f_value == 2, then *nothing is
     * done*, since the Configuration is meant to say "only save/map the dual
     * solution". In all other cases, the flow solution is saved. */

    void get_var_solution(Configuration *solc = nullptr) override
    {/*
      if (!f_Block) // no [MCF]Block to write to
        return;     // cowardly and silently return

      auto tsolc = dynamic_cast<SimpleConfiguration<int> *>(solc);
      if (tsolc && (tsolc->f_value == 2))
        return;

      auto MCFB = static_cast<MCFBlock *>(f_Block);
      MCFBlock::Vec_FNumber X(MCFB->get_NArcs());
      this->MCFGetX(X.data());
      MCFB->set_x(X.begin());
   */}
  
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    /// write the "current" dual solution in the Constraint of the MCFBlock
    /** Write the "current" dual solution, i.e., node potentials and flow
     * reduced costs, in the dual variables of the Constraint (respectively,
     * the flow conservation constraints and bound ones) of the MCFBlock. To
     * keep the same format as MCFBlock::get_Solution() and
     * MCFBlock::map[forward/back]_Modification(), the Configuration *solc can
     * be used to "partly" save it. In particular, if solc != nullptr, it is
     * a SimpleConfiguration< int >, and solc->f_value == 1, then *nothing is
     * done*, since the Configuration is meant to say "only save/map the primal
     * solution". In all other cases, the flow solution is saved. */
    
    void get_dual_solution(Configuration *solc = nullptr) override
    {
    /*  if (!f_Block) // no [MCF]Block to write to
        return;     // cowardly and silently return

      auto tsolc = dynamic_cast<SimpleConfiguration<int> *>(solc);
      if (tsolc && (tsolc->f_value == 1))
        return;

      auto MCFB = static_cast<MCFBlock *>(f_Block);
      MCFBlock::Vec_CNumber Pi(MCFB->get_NNodes());
      this->MCFGetPi(Pi.data());
      MCFB->set_pi(Pi.begin());

      MCFBlock::Vec_FNumber RC(MCFB->get_NArcs());
      this->MCFGetRC(RC.data());
      MCFB->set_rc(RC.begin());
    */}
    
    /*--------------------------------------------------------------------------*/

   // bool new_var_solution(void) override { return (this->HaveNewX()); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   // bool new_dual_solution(void) override { return (this->HaveNewPi()); }

    /*--------------------------------------------------------------------------*/
    /*
     virtual void set_unbounded_threshold( const OFValue thr ) override { }
    */

    /*--------------------------------------------------------------------------*/

    bool has_var_direction(void) override { return (true); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    bool has_dual_direction(void) override { return (true); }

    /*--------------------------------------------------------------------------*/
    /// write the current direction in the x ColVariable of the MCFBlock
    /** Write the unbounded primal direction, i.e., augmenting cycle with
     * negative cost and unbounded capacity, in the x ColVariable of the
     * MCFBlock. To keep the same format as MCFBlock::get_Solution() and
     * MCFBlock::map[forward/back]_Modification(), the Configuration *solc can
     * be used to "partly" save it. In particular, if solc != nullptr, it is
     * a SimpleConfiguration< int >, and solc->f_value == 2, then *nothing is done*,
     * since the Configuration is meant to say "only save/map the dual
     * information". In all other cases, the direction (cycle) is saved.
     *
     * Or, rather, THIS SHOULD BE DONE, BUT THE METHOD IS NOT IMPLEMENTED yet. */

    void get_var_direction(Configuration *dirc = nullptr) override
    {
      auto tsolc = dynamic_cast<SimpleConfiguration<int> *>(dirc);
      if (tsolc && (tsolc->f_value == 2))
        return;

      throw(std::logic_error(
          "MCFSolver::get_var_direction() not implemented yet"));

      // TODO: implement using MCFC::MCFGetUnbCycl()
      // anyway, unsure if any current :MCFClass properly implemente the latter
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    /// write the current dual direction in the Constraint of the MCFBlock
    /** Write the current unbounded dual direction, i.e., a cut separating two
     * shores to that the residual demand in one is greater than the capacity
     * across them, in the Constraint of the Block, in particular in the dual
     * variables of the flow conservation ones. To keep the same format as
     * MCFBlock::get_Solution() and MCFBlock::map[forward/back]_Modification(),
     * the Configuration *solc can be used to "partly" save it. In particular, if
     * solc != nullptr, it is a SimpleConfiguration< int >, and solc->f_value == 1,
     * then *nothing is done*, since the Configuration is meant to say "only
     * save/map the primal information". In all other cases, the direction (cut)
     * is saved.
     *
     * Or, rather, THIS SHOULD BE DONE, BUT THE METHOD IS NOT IMPLEMENTED yet. */

    void get_dual_direction(Configuration *dirc = nullptr) override
    {
      auto tsolc = dynamic_cast<SimpleConfiguration<int> *>(dirc);
      if (tsolc && (tsolc->f_value == 1))
        return;

      throw(std::logic_error(
          "MCFSolver::get_dual_direction() not implemented yet"));

      // TODO: implement using MCFC::MCFGetUnfCut()
      // anyway, unsure if any current :MCFClass properly implemente the latter
    }

    /*--------------------------------------------------------------------------*/
    /*
     virtual bool new_var_direction( void ) override { return( false ); }

     virtual bool new_dual_direction( void ) override{ return( false ); }
    */
    /** @} ---------------------------------------------------------------------*/
    /*-------------- METHODS FOR READING THE DATA OF THE Solver ----------------*/
    /*--------------------------------------------------------------------------*/

    /*
     virtual bool is_dual_exact( void ) const override { return( true ); }
    */

    /*--------------------------------------------------------------------------*/
    /// "publicize" MCFClass::WriteMCF
    /** Make the method
     *
     *      void WriteMCF( ostream &oStrm , int frmt = 0 )
     *
     * of the base (private) MCFClass public, so that it can be freely used. */
    //TODO: change MCFC function to Algo function.
    //using MCFC::WriteMCF;

    /*--------------------------------------------------------------------------*/
    /*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Handling the parameters of the MCFLemonSolver
     *
     * Each MCFLemonSolver< Algo > may have its own extra int / double parameters. If
     * this is the case, it will have to specialize the following methods to
     * handle them. The general definition just handles the case of the
     *
     * intLastParCDAS ==> kReopt             whether or not to reoptimize
     *
     *  @{ */


    ///@return number of int algorithmic parameters
    [[nodiscard]] idx_type get_num_int_par(void) const override
    {
        return (intLastParLEMON_NS);
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @return number of dbl algorithmic parameters
    [[nodiscard]] idx_type get_num_dbl_par(void) const override
    {
        return (dblLastParLEMON_NS);
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    /// @return number of str algorithimc parameters 
    [[nodiscard]] idx_type get_num_str_par(void) const override
    {
      return (MCFLemonSolver::strLastParLEMON);
    }
    
    /*--------------------------------------------------------------------------*/
    
    /// @brief used for obtain default value of an int parameter
    /// @param par 
    /// @return default value of parameter par, if exists
    [[nodiscard]] int get_dflt_int_par(idx_type par) const override
    {
      if( par > intLastParLEMON_NS){
        throw std::invalid_argument("Invalid int parameter: out_of_range " + std::to_string(par));
      }

      switch(par){
        case kPivot: return static_cast<int>(NetworkSimplex<GR, V, C>::PivotRule::BLOCK_SEARCH);
        default: return(CDASolver::get_dflt_int_par(par));
        
      }
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for obtain default value of a double parameter
    /// @param par 
    /// @return default value of parameter par, if exists
    [[nodiscard]] double get_dflt_dbl_par(idx_type par) const override
    {
      if( par > dblLastParLEMON_NS){
        throw std::invalid_argument("Invalid dbl parameter: out_of_range " + std::to_string(par));
      }
      switch(par){
        default: return(CDASolver::get_dflt_dbl_par(par)); 
      }

      
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    /// @brief used for obtain default value of a string parameter
    /// @param par
    /// @return default valule of parameter par, if exists
    [[nodiscard]] const std::string &get_dflt_str_par(idx_type par)
        const override
    {
      if(par > MCFLemonSolverBase<NetworkSimplex, GR, V, C>::strLastParLEMON){
        throw std::invalid_argument("Invalid str parameter: out_of_range " + std::to_string(par));
      }
      static const std::string _empty;
      if (par == MCFLemonSolverBase<NetworkSimplex, GR, V, C>::strLastParLEMON)
        return (_empty);

      return (CDASolver::get_dflt_str_par(par));
    }
        

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    /// @brief used for get the value of string parameters
    /// @param par 
    /// @return value of parameter indexed by par
    [[nodiscard]] const std::string & get_str_par( idx_type par ) const override {
        
        if( par == MCFLemonSolverBase<NetworkSimplex, GR, V, C>::strDMXFile )
        return( this->f_dmx_file );

        return( get_dflt_str_par( par ) );
        }
   

    /*--------------------------------------------------------------------------*/
    
    /// @brief used for get the value of int parameters
    /// @param par 
    /// @return value of parameter indexed by par
    [[nodiscard]] int get_int_par(idx_type par) const override
    {
    
      if(par == kPivot){
        return f_pivot_rule_type;
      }

      return (get_dflt_int_par(par));
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for get the value of dbl parameters
    /// @param par 
    /// @return value of parameter indexed by par
    [[nodiscard]] double get_dbl_par(idx_type par) const override
    {
      //Da finire parametri algoritmici dbl

      return (get_dflt_dbl_par(par));
    }
        
    /*--------------------------------------------------------------------------*/
    
    /// @brief used for convert string to int parameter's index
    /// @param name 
    /// @return an index that denotes parameter name
    [[nodiscard]] idx_type int_par_str2idx(const std::string &name)
        const override
    {
      if (name == "kPivot")
        return (kPivot);

      return (CDASolver::int_par_str2idx(name));
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for convert string to dbl parameter's index
    /// @param name 
    /// @return an index that denotes parameter name
    [[nodiscard]] idx_type dbl_par_str2idx(const std::string &name)
        const override
    {
      return (CDASolver::dbl_par_str2idx(name));
    }
      
    /*--------------------------------------------------------------------------*/
    
    /// @brief used for convert int index idx to phrasal rapresentation
    /// @param idx 
    /// @return a string that denotes index idx parameter
    [[nodiscard]] const std::string &int_par_idx2str(idx_type idx)
        const override
    {

      if(idx > intLastParLEMON_NS){
        throw std::invalid_argument("Invalid index: out_of_range " + std::to_string(idx));
      }
      static const std::string par = "kPivot";
      switch(idx){
        case kPivot: return (par);
        default: break;
      }

      return (CDASolver::int_par_idx2str(idx));
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for convert dbl index idx to phrasal rapresentation
    /// @param idx 
    /// @return a string that denotes index idx parameter
    [[nodiscard]] const std::string &dbl_par_idx2str(idx_type idx)
        const override
    {
      if(idx > dblLastParLEMON_NS){
        throw std::invalid_argument("Invalid index: out_of_range " + std::to_string(idx));
      }

      switch(idx){
        default: break;
      }

      return (CDASolver::dbl_par_idx2str(idx));
    }
    
    /** @} ---------------------------------------------------------------------*/
    /*------------ METHODS FOR HANDLING THE State OF THE MCFLemonSolver -------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Handling the State of the MCFSolver
     *  @{ */

    //[[nodiscard]] State *get_State(void) const override;

    /*--------------------------------------------------------------------------*/

    //void put_State(const State &state) override;

    /*--------------------------------------------------------------------------*/

    //void put_State(State &&state) override;

    /** @} ---------------------------------------------------------------------*/
    /*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Changing the data of the model
     *  @{ */

    /** The only reason why MCFSolver::add_Modification() needs be defined is to
     * properly react to NBModification. Indeed, the correct reaction is to
     * *immediately* reload the MCFBlock, besides clearing the list of
     * Modification as Solver::add_Modification() already does. The issue is
     * that if arcs/nodes are added/deleted after the NBModification is issued
     * but before it is processed, then the number of nodes/arcs at the moment
     * in which the NBModification is processed is different from that at the
     * moment in which is issued, which may break the "naming convention"
     * (because the name of, say, a newly created arc depends on the current
     * state and/or number of the arcs).
     *
     * Important note: THIS VERSION ONLY WORKS PROPERLY IF THE MCFBlock IS
     * "FRESHLY MINTED", I.E., THERE ARE NO CLOSED OR DELETED ARCS.
     *
     * This should ordinarily always happen, as whenever the MCFBlock is changed
     * the NBModification is immediately issued. The problem may come if the
     * MCFBlock is a R3Block of another MCFBlock which is loaded and then
     * further modified, and the NBModification to this MCFBlock is generated by
     * a map_forward_Modification() of the NBModification to the original
     * MCFBlock: then, this MCFBlock may be copied from a MCFBlock that has
     * closed or deleted arcs and this method would not work. */

    /*void add_Modification(sp_Mod &mod) override
    {
      if (std::dynamic_pointer_cast<const NBModification>(mod))
      {
        // this is the "nuclear option": the MCFBlock has been re-loaded, so
        // the MCFClass solver also has to (immediately)
        //TODO: change MCFC function to Algo function.
        auto MCFB = static_cast<MCFBlock *>(f_Block);
        MCFC::LoadNet(MCFB->get_MaxNNodes(), MCFB->get_MaxNArcs(),
                      MCFB->get_NNodes(), MCFB->get_NArcs(),
                      MCFB->get_U().empty() ? nullptr : MCFB->get_U().data(),
                      MCFB->get_C().empty() ? nullptr : MCFB->get_C().data(),
                      MCFB->get_B().empty() ? nullptr : MCFB->get_B().data(),
                      MCFB->get_SN().data(), MCFB->get_EN().data());
        // TODO: PreProcess() changes the internal data of the MCFSolver using
        //       information about how the data of the MCF is *now*. If the
        //       data changes, some of the deductions (say, reducing the capacity
        //       of some arcs) may no longer be correct and they should be undone,
        //       but there isn't any proper way to handle this. Thus, PreProcess()
        //       has to be disabled for now; maybe later on someone will take
        //       care to make this work (or maybe not).
        // MCFC::PreProcess();
        // besides, any outstanding modification makes no sense any longer
        mod_clear();
      }
      else
        push_back(mod);
    }*/

    /*--------------------------------------------------------------------------*/
    /*-------------------------------- FRIENDS ---------------------------------*/
    /*--------------------------------------------------------------------------*/

    friend class MCFLemonState; // make MCFSolverState friend

    /** @} ---------------------------------------------------------------------*/
    /*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
    /*--------------------------------------------------------------------------*/

  protected:
    /*--------------------------------------------------------------------------*/
    /*-------------------------- PROTECTED METHODS -----------------------------*/
    /*--------------------------------------------------------------------------*/

    void process_outstanding_Modification(void);

    void guts_of_poM(c_p_Mod mod);

    /*--------------------------------------------------------------------------*/
    /*---------------------------- PROTECTED FIELDS  ---------------------------*/
    /*--------------------------------------------------------------------------*/

    const static std::vector<int> Solver_2_MCFClass_int;
    // the (static const) map between Solver int parameters and MCFClass ones

    const static std::vector<int> Solver_2_MCFClass_dbl;
    // the (static const) map between Solver int parameters and MCFClass ones

    std::string f_dmx_file; 
    // string for DMX file output


    /*--------------------------------------------------------------------------*/
    /*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
    /*--------------------------------------------------------------------------*/

  private:
    /*--------------------------------------------------------------------------*/
    /*-------------------------- PRIVATE METHODS -------------------------------*/
    /*--------------------------------------------------------------------------*/

    SMSpp_insert_in_factory_h;

     /*--------------------------------------------------------------------------*/
    /*-------------------------- PRIVATE FIELDS -------------------------------*/
    /*--------------------------------------------------------------------------*/
    int status = SMSpp_di_unipi_it::LEMON_sol_type::UNSOLVED;  //Variable used in compute function for getting status
    typename NetworkSimplex< GR , V , C >::ProblemType status_2_pType;
    //Status of compute() method
    
    int f_pivot_rule_type; // Variable used in set_par for switching value with f_pivot_rule

    double ticks;  //Elaped time in ticks for compute() method
    /*--------------------------------------------------------------------------*/

  }; // end( class MCFLemonSolver<NetworkSimplex, GR, V, C>  Specialization)


/*--------------------------------------------------------------------------*/
/*-------------------------------- SPECIALIZED CLASSES ---------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------------- CYCLECANCELING --------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/


/** Specialized MCFLemonSolver<CycleCanceling, GR, V, C> that contains specialized
 *  compute() method, enums for indexing algorithimc parameters and function
 *  set/get_*_par for manage them.
 *  This class derives from MCFLemonSolverBase for using set_Block() methods and
 *  for inherit f_algo and dgp fields.
 * 
 *  Template parameters are:
 * 
 *  -  CycleCanceling implements three different cycle-canceling algorithms for finding 
 *     a minimum cost flow. The most efficent one is the Cancel-and-tighten algorithm,
 *     thus it is the default method. It runs in strongly polynomial time, but in practice, 
 *     it is typically orders of magnitude slower than the scaling algorithms and NetworkSimplex.
 * 
 *  - GR that represents the directed graph, the possibilities are described in the
 *    file MCFLemonSolver.h at line 339
 *  
 *  - V, which is the type of flows / deficits; typically, double can be used
 *    for maximum compatibility, but int (or even smaller) would yeld better
 *    performances;
 *  
 * - C, which is the type of ar costs; typically, double can be used for
 *    maximum compatibility, but int (or even smaller) would yeld better
 *    performances;
*/
  template< typename GR, typename V, typename C>
class MCFLemonSolver<CycleCanceling, GR, V, C> : virtual public CDASolver, public MCFLemonSolverBase<CycleCanceling, GR, V, C>, Fields< CycleCanceling<GR, V, C> >
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

   

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
   kUnEval = 0     compute() has not been called yet

   kUnbounded = kUnEval + 1     the model is provably unbounded
 *  @{ */
  const int kErrorStatus = -1;
/*
    
    kUnEval = 0     compute() has not been called yet

    kUnbounded = kUnEval + 1     the model is provably unbounded

    kInfeasible                  the model is provably infeasible

    kBothInfeasible = kInfeasible + 1     both primal and dual infeasible

    kOK = 7         successful compute()
                    Any return value between kUnEval (excluded) and kOK
        (included) means that the object ran smoothly

    kStopTime = kOK + 1          stopped because of time limit

    kStopIter                    stopped because of iteration limit

    kError = 15     compute() stopped because of unrecoverable error
                    Any return value >= kError means that the object was
         forced to stop due to some error, e.g. of numerical nature

    kLowPrecision = kError + 1   a solution found but not provably optimal
    */



/** @} ---------------------------------------------------------------------*/
/*-------------- CONSTRUCTING AND DESTRUCTING MCFLemonSolver ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing MCFLemonSolver
 *  @{ */

 /// constructor: assign the default parameter
 /** Void constructor: Build f_method_type with default parameter */

 MCFLemonSolver( void ) : CDASolver(), MCFLemonSolverBase<CycleCanceling, GR, V, C>(),  Fields<CycleCanceling<GR, V, C > >() {
  f_method_type = CycleCanceling<GR, V, C>::Method::CANCEL_AND_TIGHTEN;
  }
 
 ///destructor:
 /**Void destructor: Delete f_method from memory */
  ~MCFLemonSolver( void ) {
    delete Fields<CycleCanceling< GR, V, C > >::f_method;   
  }
 
 

/*--------------------------------------------------------------------------*/
 /* intMaxIter = 0     maximum iterations for the next call to solve()

    intMaxSol          maximum number of different solutions to report

    intLogVerb         "verbosity" of the log

    intMaxDSol         maximum number of different dual solutions

    intLastParCDAS     first allowed parameter value for derived classes
    */

 enum LEMON_CC_int_par_type{
        kMethod = intLastParCDAS,
        intLastParLEMON_CC ///< first allowed parameter value for derived classes
                           /**< convenience value for easily allow derived classes
                            * to further extend the set of types of return codes */
 };

/*--------------------------------------------------------------------------*/
 /* dblMaxTime = 0    maximum time for the next call to solve()

    dblRelAcc         relative accuracy for declaring a solution optimal

    dblAbsAcc          absolute accuracy for declaring a solution optimal

    dblUpCutOff        upper cutoff for stopping the algorithm

    dblLwCutOff        lower cutoff for stopping the algorithm

    dblRAccSol          maximum relative error in any reported solution

    dblAAccSol          maximum absolute error in any reported solution

    dblFAccSol          maximum constraint violation in any reported solution

    dblRAccDSol         maximum relative error in any dual solution

    dblAAccDSol         maximum absolute error in any dual solution

    dblFAccDSol         maximum absolute error in any dual solution

    dblLastParCDAS      first allowed parameter value for derived classes
    */

/*--------------------------------------------------------------------------*/

/// public enum for the type of the solution
/// to MCFLemonSolver< CycleCanceling , C , V >

enum dbl_par_type_LEMON_CC{
  dblCycleCancelingFactor = dblLastParCDAS ,
  dblLastParLEMON_CC ///< first allowed parameter value for derived classes
                           /**< convenience value for easily allow derived classes
                            * to further extend the set of types of return codes */
};

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *
 * Parameter-wise, MCFSolver maps the parameters of [CDA]Solver
 *
 *  intMaxIter = 0    maximum iterations for the next call to solve()
 *  intMaxSol         maximum number of different solutions to report
 *  intLogVerb        "verbosity" of the log
 *  intMaxDSol        maximum number of different dual solutions
 *
 *  dblMaxTime = 0    maximum time for the next call to solve()
 *  dblRelAcc         relative accuracy for declaring a solution optimal
 *  dblAbsAcc         absolute accuracy for declaring a solution optimal
 *  dblUpCutOff       upper cutoff for stopping the algorithm
 *  dblLwCutOff       lower cutoff for stopping the algorithm
 *  dblRAccSol        maximum relative error in any reported solution
 *  dblAAccSol        maximum absolute error in any reported solution
 *  dblFAccSol        maximum constraint violation in any reported solution
 *  dblRAccDSol       maximum relative error in any dual solution
 *  dblAAccDSol       maximum absolute error in any dual solution
 *  dblFAccDSol       maximum absolute error in any dual solution
 *
 * into the parameter of MCFClass
 *
 * kMaxTime = 0       max time
 * kMaxIter           max number of iteration
 * kEpsFlw            tolerance for flows
 * kEpsDfct           tolerance for deficits
 * kEpsCst            tolerance for costs
 *
 * It then "extends" them, using
 *
 *  intLastParCDAS    first allowed parameter value for derived classes
 *  dblLastParCDAS    first allowed parameter value for derived classes
 *
 * In particular, one now has
 *
 * intLastParCDAS ==> kReopt             whether or not to reoptimize
 *
 * and any other parameter of specific :MCFClass following. This is done
 * via the two const static arrays Solver_2_MCFClass_int and
 * Solver_2_MCFClass_dbl, with a negative entry meaning "there is no such
 * parameter in MCFSolver".
 *
 *  @{ */

 

    /*--------------------------------------------------------------------------*/
    // set the ostream for the Solver log
    // not really, MCFClass objects are remarkably silent
    //
    // virtual void set_log( std::ostream *log_stream = nullptr ) override;

    /*--------------------------------------------------------------------------*/
    /// @brief set the parameter par with value
    /// @param par 
    /// @param value 
    void set_par(idx_type par, int value){
        if(par != kMethod){
                CDASolver::set_par(par, value);
                return;
        }

        if((value < 0) || (value > 4)){
                throw(std::invalid_argument(std::to_string(value)));
        }

        if(value == f_method_type){
                return; //Nothing is changed
        }

        f_method_type = value;
        delete Fields<CycleCanceling<GR, V, C>>::f_method;
        Fields<CycleCanceling<GR, V, C>>::f_method = NULL;
        return;
    } 
 
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    //TODO: change MCFC function to Algo function.
    /*void set_par(idx_type par, double value) override
    {
      if (Solver_2_MCFClass_dbl[par] >= 0)
      ;
      //  MCFC::SetPar(Solver_2_MCFClass_dbl[par], double(value));
    }*/
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/


    /** @} ---------------------------------------------------------------------*/
    /*--------------------- METHODS FOR SOLVING THE Block ----------------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Solving the MCF encoded by the current MCFBlock
     *  @{ */
    /// (try to) solve the MCF encoded in the MCFBlock 
        int compute( bool changedvars = true ) override {
        const static std::array<int, 6> LemonStatus_2_MCFstatus = {
                kErrorStatus, SMSpp_di_unipi_it::LEMON_sol_type::OPTIMAL,
                 kErrorStatus , SMSpp_di_unipi_it::LEMON_sol_type::INFEASIBLE,
                SMSpp_di_unipi_it::LEMON_sol_type::UNBOUNDED, kErrorStatus};
        
        const static std::array<int, 6> MCFstatus_2_sol_type = {
                kUnEval, Solver::kOK, kStopTime, kInfeasible, Solver::kUnbounded,
                Solver::kError};

        lock(); // first of all, acquire self-lock

        if (!f_Block)            // there is no [MCFBlock] to solve
                return (kBlockLocked); // return error

        bool owned = f_Block->is_owned_by(f_id); // check if already locked
        if ((!owned) && (!f_Block->read_lock())) // if not try to read_lock
                return (kBlockLocked);                 // return error on failure

        // while [read_]locked, process any outstanding Modification
        //TODO: ensure that modification are actually processed for MCFLemonSolver.
        //process_outstanding_Modification();

        if (!f_dmx_file.empty())
        { // if so required
                // output the current instance (after the changes) to a DMX file
                std::ofstream ProbFile(f_dmx_file, ios_base::out | ios_base::trunc);
                if (!ProbFile.is_open())
                throw(std::logic_error("cannot open DMX file " + f_dmx_file));

                //WriteMCF(ProbFile);
                writeDimacsMat(ProbFile, *MCFLemonSolver::dgp);
                ProbFile.close();
        }

        if (!owned)               // if the [MCF]Block was actually read_locked
                f_Block->read_unlock(); // read_unlock it

        auto start = chrono::system_clock::now();
        
        if(Fields<CycleCanceling< GR, V, C > >::f_method != NULL){
        this->status = MCFLemonSolver::f_algo->run(*Fields<CycleCanceling< GR, V, C > >::f_method);
        }else{
                //Build f_method and execute run() method
                Fields<CycleCanceling<GR, V, C> >::f_method = new typename   CycleCanceling<GR, V, C>::Method(f_method_type);
                this->status = MCFLemonSolver::f_algo->run(*Fields<CycleCanceling< GR, V, C > >::f_method);
               
        }
        
        auto end = chrono::system_clock::now();

        chrono::duration< double > elapsed = end - start;
        ticks = elapsed.count();


        
        
        unlock(); // release self-lock

        // now give out the result: note that the vector MCFstatus_2_sol_type[]
        // starts from 0 whereas the first value of MCFStatus is -1 (= kUnSolved),
        // hence the returned status has to be shifted by + 1
        return (MCFstatus_2_sol_type[this->get_status()]);
        }

    

    /** @} ---------------------------------------------------------------------*/
    /*---------------------- METHODS FOR READING RESULTS -----------------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Accessing the found solutions (if any)
     *  @{ */
    
    

    ///Return the status of the run() in compute method
    int get_status(void) const 
    {
      return (this->status);
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    //Return the elapsed time for the run() in compute method   
    double get_elapsed_time(void) const override
    {
      return (this->ticks);
    }



    /*--------------------------------------------------------------------------*/
    //Return the lower bound solution(optimal) for the problem
    OFValue get_lb(void) override { return OFValue(MCFLemonSolver::f_algo->totalCost()); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    //Return the upper bound solution(optimal) for the problem
    OFValue get_ub(void) override { return OFValue(MCFLemonSolver::f_algo->totalCost()); }

    /*--------------------------------------------------------------------------*/
    
    /// @brief Analyze the return status of run() 
    /// @param  
    /// @return if the solution found is OPTIMAL or UNBOUNDED returns true else false
    bool has_var_solution(void) override
    {
     switch (this->get_status()) {
      case( CycleCanceling< GR , V , C >::ProblemType::OPTIMAL ):
      case( CycleCanceling< GR , V , C >::ProblemType::UNBOUNDED ):  
       return( true );
      default:
       return( false );
      }
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief Analyze the return status of run()
    /// @param  
    /// @return if the solution found is OPTIMAL or INFEASIBLE returns true else false
    bool has_dual_solution(void) override
    {
     switch( this->get_status() ) {
      case( CycleCanceling< GR , V , C >::ProblemType::OPTIMAL ):
      case( CycleCanceling< GR , V , C >::ProblemType::INFEASIBLE ):
       return( true );
      default:
       return( false );
      }
    }

    /*--------------------------------------------------------------------------*/
    /*
     virtual bool is_var_feasible( void ) override { return( true ); }

     virtual bool is_dual_feasible( void ) override { return( true ); }
    */
    /*--------------------------------------------------------------------------*/
    /// write the "current" flow in the x ColVariable of the MCFBlock
    /** Write the "current" flow in the x ColVariable of the MCFBlock. To keep
     * the same format as MCFBlock::get_Solution() and
     * MCFBlock::map[forward/back]_Modification(), the Configuration *solc can
     * be used to "partly" save it. In particular, if solc != nullptr, it is
     * a SimpleConfiguration< int >, and solc->f_value == 2, then *nothing is
     * done*, since the Configuration is meant to say "only save/map the dual
     * solution". In all other cases, the flow solution is saved. */

    void get_var_solution(Configuration *solc = nullptr) override
    {/*
      if (!f_Block) // no [MCF]Block to write to
        return;     // cowardly and silently return

      auto tsolc = dynamic_cast<SimpleConfiguration<int> *>(solc);
      if (tsolc && (tsolc->f_value == 2))
        return;

      auto MCFB = static_cast<MCFBlock *>(f_Block);
      MCFBlock::Vec_FNumber X(MCFB->get_NArcs());
      this->MCFGetX(X.data());
      MCFB->set_x(X.begin());
   */}
  
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    /// write the "current" dual solution in the Constraint of the MCFBlock
    /** Write the "current" dual solution, i.e., node potentials and flow
     * reduced costs, in the dual variables of the Constraint (respectively,
     * the flow conservation constraints and bound ones) of the MCFBlock. To
     * keep the same format as MCFBlock::get_Solution() and
     * MCFBlock::map[forward/back]_Modification(), the Configuration *solc can
     * be used to "partly" save it. In particular, if solc != nullptr, it is
     * a SimpleConfiguration< int >, and solc->f_value == 1, then *nothing is
     * done*, since the Configuration is meant to say "only save/map the primal
     * solution". In all other cases, the flow solution is saved. */
    
    void get_dual_solution(Configuration *solc = nullptr) override
    {
    /*  if (!f_Block) // no [MCF]Block to write to
        return;     // cowardly and silently return

      auto tsolc = dynamic_cast<SimpleConfiguration<int> *>(solc);
      if (tsolc && (tsolc->f_value == 1))
        return;

      auto MCFB = static_cast<MCFBlock *>(f_Block);
      MCFBlock::Vec_CNumber Pi(MCFB->get_NNodes());
      this->MCFGetPi(Pi.data());
      MCFB->set_pi(Pi.begin());

      MCFBlock::Vec_FNumber RC(MCFB->get_NArcs());
      this->MCFGetRC(RC.data());
      MCFB->set_rc(RC.begin());
    */}
    
    /*--------------------------------------------------------------------------*/

   // bool new_var_solution(void) override { return (this->HaveNewX()); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   // bool new_dual_solution(void) override { return (this->HaveNewPi()); }

    /*--------------------------------------------------------------------------*/
    /*
     virtual void set_unbounded_threshold( const OFValue thr ) override { }
    */

    /*--------------------------------------------------------------------------*/

    bool has_var_direction(void) override { return (true); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    bool has_dual_direction(void) override { return (true); }

    /*--------------------------------------------------------------------------*/
    /// write the current direction in the x ColVariable of the MCFBlock
    /** Write the unbounded primal direction, i.e., augmenting cycle with
     * negative cost and unbounded capacity, in the x ColVariable of the
     * MCFBlock. To keep the same format as MCFBlock::get_Solution() and
     * MCFBlock::map[forward/back]_Modification(), the Configuration *solc can
     * be used to "partly" save it. In particular, if solc != nullptr, it is
     * a SimpleConfiguration< int >, and solc->f_value == 2, then *nothing is done*,
     * since the Configuration is meant to say "only save/map the dual
     * information". In all other cases, the direction (cycle) is saved.
     *
     * Or, rather, THIS SHOULD BE DONE, BUT THE METHOD IS NOT IMPLEMENTED yet. */

    void get_var_direction(Configuration *dirc = nullptr) override
    {
      auto tsolc = dynamic_cast<SimpleConfiguration<int> *>(dirc);
      if (tsolc && (tsolc->f_value == 2))
        return;

      throw(std::logic_error(
          "MCFSolver::get_var_direction() not implemented yet"));

      // TODO: implement using MCFC::MCFGetUnbCycl()
      // anyway, unsure if any current :MCFClass properly implemente the latter
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    /// write the current dual direction in the Constraint of the MCFBlock
    /** Write the current unbounded dual direction, i.e., a cut separating two
     * shores to that the residual demand in one is greater than the capacity
     * across them, in the Constraint of the Block, in particular in the dual
     * variables of the flow conservation ones. To keep the same format as
     * MCFBlock::get_Solution() and MCFBlock::map[forward/back]_Modification(),
     * the Configuration *solc can be used to "partly" save it. In particular, if
     * solc != nullptr, it is a SimpleConfiguration< int >, and solc->f_value == 1,
     * then *nothing is done*, since the Configuration is meant to say "only
     * save/map the primal information". In all other cases, the direction (cut)
     * is saved.
     *
     * Or, rather, THIS SHOULD BE DONE, BUT THE METHOD IS NOT IMPLEMENTED yet. */

    void get_dual_direction(Configuration *dirc = nullptr) override
    {
      auto tsolc = dynamic_cast<SimpleConfiguration<int> *>(dirc);
      if (tsolc && (tsolc->f_value == 1))
        return;

      throw(std::logic_error(
          "MCFSolver::get_dual_direction() not implemented yet"));

      // TODO: implement using MCFC::MCFGetUnfCut()
      // anyway, unsure if any current :MCFClass properly implemente the latter
    }

    /*--------------------------------------------------------------------------*/
    /*
     virtual bool new_var_direction( void ) override { return( false ); }

     virtual bool new_dual_direction( void ) override{ return( false ); }
    */
    /** @} ---------------------------------------------------------------------*/
    /*-------------- METHODS FOR READING THE DATA OF THE Solver ----------------*/
    /*--------------------------------------------------------------------------*/

    /*
     virtual bool is_dual_exact( void ) const override { return( true ); }
    */

    /*--------------------------------------------------------------------------*/
    /// "publicize" MCFClass::WriteMCF
    /** Make the method
     *
     *      void WriteMCF( ostream &oStrm , int frmt = 0 )
     *
     * of the base (private) MCFClass public, so that it can be freely used. */
    //TODO: change MCFC function to Algo function.
    //using MCFC::WriteMCF;

    /*--------------------------------------------------------------------------*/
    /*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Handling the parameters of the MCFLemonSolver
     *
     * Each MCFLemonSolver< Algo > may have its own extra int / double parameters. If
     * this is the case, it will have to specialize the following methods to
     * handle them. The general definition just handles the case of the
     *
     * intLastParCDAS ==> kReopt             whether or not to reoptimize
     *
     *  @{ */


     ///@return number of int algorithmic parameters
    [[nodiscard]] idx_type get_num_int_par(void) const override
    {
        return (intLastParLEMON_CC);
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @return number of dbl algorithmic parameters
    [[nodiscard]] idx_type get_num_dbl_par(void) const override
    {
        return (dblLastParLEMON_CC);
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @return number of str algorithimc parameters 
    [[nodiscard]] idx_type get_num_str_par(void) const override
    {
      return (MCFLemonSolver::strLastParLEMON);
    }
    
    /*--------------------------------------------------------------------------*/
    
    /// @brief used for obtain default value of an int parameter
    /// @param par 
    /// @return default value of parameter par, if exists
    [[nodiscard]] int get_dflt_int_par(idx_type par) const override
    {
      if(par > intLastParLEMON_CC){
        throw(std::invalid_argument(std::to_string(par)));
      }

      switch(par){
        case kMethod: return CycleCanceling<GR, V, C>::Method::CANCEL_AND_TIGHTEN;
        default: return (CDASolver::get_dflt_int_par(par));
      }
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for obtain default value of a double parameter
    /// @param par 
    /// @return default value of parameter par, if exists
    [[nodiscard]] double get_dflt_dbl_par(idx_type par) const override
    {
      if(par > dblLastParLEMON_CC){
        throw(std::invalid_argument(std::to_string(par)));
      }

      switch(par){
        default: return (CDASolver::get_dflt_dbl_par(par));
      }
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for obtain default value of a string parameter
    /// @param par
    /// @return default valule of parameter par, if exists
    [[nodiscard]] const std::string &get_dflt_str_par(idx_type par)
        const override
    {

       if(par > MCFLemonSolverBase<CycleCanceling, GR, V, C>::strLastParLEMON){
        throw std::invalid_argument("Invalid str parameter: out_of_range " + std::to_string(par));
      }

      static const std::string _empty;
      if (par == MCFLemonSolverBase<CycleCanceling, GR, V, C>::strLastParLEMON)
        return (_empty);

      return (CDASolver::get_dflt_str_par(par));
    }
    
    /*--------------------------------------------------------------------------*/
    
    /// @brief used for get the value of int parameters
    /// @param par 
    /// @return value of parameter indexed by par
    [[nodiscard]] int get_int_par(idx_type par) const override
    {
      
      if(par == kMethod){
        return f_method_type;
      }

      return (get_dflt_int_par(par));
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for get the value of dbl parameters
    /// @param par 
    /// @return value of parameter indexed by par
    [[nodiscard]] double get_dbl_par(idx_type par) const override
    {
      //Da finire parametri algoritmici dbl
      return (get_dflt_dbl_par(par));
    }
    
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    /// @brief used for get the value of string parameters
    /// @param par 
    /// @return value of parameter indexed by par
    [[nodiscard]] const std::string & get_str_par( idx_type par ) const override {
        
      if( par == MCFLemonSolverBase<CycleCanceling, GR, V, C>::strDMXFile )
        return( this->f_dmx_file );

      return( get_dflt_str_par( par ) );
    }
   

    /*--------------------------------------------------------------------------*/
    
    /// @brief used for convert string to int parameter's index
    /// @param name 
    /// @return an index that denotes parameter name
    [[nodiscard]] idx_type int_par_str2idx(const std::string &name)
        const override
    {
      if (name == "kMethod")
        return (kMethod);

      return (CDASolver::int_par_str2idx(name));
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for convert string to dbl parameter's index
    /// @param name 
    /// @return an index that denotes parameter name
    [[nodiscard]] idx_type dbl_par_str2idx(const std::string &name)
        const override
    {
      return (CDASolver::dbl_par_str2idx(name));
    }
      
    /*--------------------------------------------------------------------------*/
    
    /// @brief used for convert int index idx to phrasal rapresentation
    /// @param idx 
    /// @return a string that denotes index idx parameter
    [[nodiscard]] const std::string &int_par_idx2str(idx_type idx)
        const override
    {

      if(idx > intLastParLEMON_CC){
        throw std::invalid_argument(std::to_string(idx));
      }
      static const std::string par = "kMethod";
      switch(idx){
        case kMethod: return (par);
        default: break;
      }


      return (CDASolver::int_par_idx2str(idx));
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for convert dbl index idx to phrasal rapresentation
    /// @param idx 
    /// @return a string that denotes index idx parameter
    [[nodiscard]] const std::string &dbl_par_idx2str(idx_type idx)
        const override
    {
      if(idx > dblLastParLEMON_CC){
        throw std::invalid_argument(std::to_string(idx));
      }

      switch(idx){
        default: break;
      }

      return (CDASolver::dbl_par_idx2str(idx));
    }
    
    /** @} ---------------------------------------------------------------------*/
    /*------------ METHODS FOR HANDLING THE State OF THE MCFLemonSolver -------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Handling the State of the MCFSolver
     *  @{ */

    //[[nodiscard]] State *get_State(void) const override;

    /*--------------------------------------------------------------------------*/

    //void put_State(const State &state) override;

    /*--------------------------------------------------------------------------*/

    //void put_State(State &&state) override;

    /** @} ---------------------------------------------------------------------*/
    /*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Changing the data of the model
     *  @{ */

    /** The only reason why MCFSolver::add_Modification() needs be defined is to
     * properly react to NBModification. Indeed, the correct reaction is to
     * *immediately* reload the MCFBlock, besides clearing the list of
     * Modification as Solver::add_Modification() already does. The issue is
     * that if arcs/nodes are added/deleted after the NBModification is issued
     * but before it is processed, then the number of nodes/arcs at the moment
     * in which the NBModification is processed is different from that at the
     * moment in which is issued, which may break the "naming convention"
     * (because the name of, say, a newly created arc depends on the current
     * state and/or number of the arcs).
     *
     * Important note: THIS VERSION ONLY WORKS PROPERLY IF THE MCFBlock IS
     * "FRESHLY MINTED", I.E., THERE ARE NO CLOSED OR DELETED ARCS.
     *
     * This should ordinarily always happen, as whenever the MCFBlock is changed
     * the NBModification is immediately issued. The problem may come if the
     * MCFBlock is a R3Block of another MCFBlock which is loaded and then
     * further modified, and the NBModification to this MCFBlock is generated by
     * a map_forward_Modification() of the NBModification to the original
     * MCFBlock: then, this MCFBlock may be copied from a MCFBlock that has
     * closed or deleted arcs and this method would not work. */

    /*void add_Modification(sp_Mod &mod) override
    {
      if (std::dynamic_pointer_cast<const NBModification>(mod))
      {
        // this is the "nuclear option": the MCFBlock has been re-loaded, so
        // the MCFClass solver also has to (immediately)
        //TODO: change MCFC function to Algo function.
        auto MCFB = static_cast<MCFBlock *>(f_Block);
        MCFC::LoadNet(MCFB->get_MaxNNodes(), MCFB->get_MaxNArcs(),
                      MCFB->get_NNodes(), MCFB->get_NArcs(),
                      MCFB->get_U().empty() ? nullptr : MCFB->get_U().data(),
                      MCFB->get_C().empty() ? nullptr : MCFB->get_C().data(),
                      MCFB->get_B().empty() ? nullptr : MCFB->get_B().data(),
                      MCFB->get_SN().data(), MCFB->get_EN().data());
        // TODO: PreProcess() changes the internal data of the MCFSolver using
        //       information about how the data of the MCF is *now*. If the
        //       data changes, some of the deductions (say, reducing the capacity
        //       of some arcs) may no longer be correct and they should be undone,
        //       but there isn't any proper way to handle this. Thus, PreProcess()
        //       has to be disabled for now; maybe later on someone will take
        //       care to make this work (or maybe not).
        // MCFC::PreProcess();
        // besides, any outstanding modification makes no sense any longer
        mod_clear();
      }
      else
        push_back(mod);
    }*/

    /*--------------------------------------------------------------------------*/
    /*-------------------------------- FRIENDS ---------------------------------*/
    /*--------------------------------------------------------------------------*/

    friend class MCFLemonState; // make MCFSolverState friend

    /** @} ---------------------------------------------------------------------*/
    /*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
    /*--------------------------------------------------------------------------*/

  protected:
    /*--------------------------------------------------------------------------*/
    /*-------------------------- PROTECTED METHODS -----------------------------*/
    /*--------------------------------------------------------------------------*/

    void process_outstanding_Modification(void);

    void guts_of_poM(c_p_Mod mod);

    /*--------------------------------------------------------------------------*/
    /*---------------------------- PROTECTED FIELDS  ---------------------------*/
    /*--------------------------------------------------------------------------*/

    const static std::vector<int> Solver_2_MCFClass_int;
    // the (static const) map between Solver int parameters and MCFClass ones

    const static std::vector<int> Solver_2_MCFClass_dbl;
    // the (static const) map between Solver int parameters and MCFClass ones

    std::string f_dmx_file; 
    // string for DMX file output


    /*--------------------------------------------------------------------------*/
    /*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
    /*--------------------------------------------------------------------------*/

  private:
    /*--------------------------------------------------------------------------*/
    /*-------------------------- PRIVATE METHODS -------------------------------*/
    /*--------------------------------------------------------------------------*/

    SMSpp_insert_in_factory_h;

     /*--------------------------------------------------------------------------*/
    /*-------------------------- PRIVATE FIELDS -------------------------------*/
    /*--------------------------------------------------------------------------*/

    int status = SMSpp_di_unipi_it::LEMON_sol_type::UNSOLVED;  //Variable used in compute function for getting status
    typename CycleCanceling< GR , V , C >::ProblemType status_2_pType;
    //Status of compute() method
    int f_method_type; //Variable used in set_par for switching value with f_method
    double ticks;  //Elaped time in ticks for compute() method
    /*--------------------------------------------------------------------------*/

  }; // end( class MCFLemonSolver<CycleCanceling> specialization )

/*--------------------------------------------------------------------------*/
/*-------------------------------- SPECIALIZED CLASSES ---------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------------- CAPACITYSCALING -------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/


/** Specialized MCFLemonSolver<CapacityScaling, GR, V, C> that contains specialized
 *  compute() method, enums for indexing algorithimc parameters and function
 *  set/get_*_par for manage them.
 *  This class derives from MCFLemonSolverBase for using set_Block() methods and
 *  for inherit f_algo and dgp fields.
 * 
 *  Template parameters are:
 * 
 *   - CapacityScaling implements the capacity scaling version of the successive shortest path
 *     algorithm for finding a minimum cost flow. It is an efficient dual solution method,
 *     which runs in polynomial time.
 *     In special case it can be more efficient than CostScaling and NetworkSimplex algorithms.
 * 
 *  - GR that represents the directed graph, the possibilities are described in the
 *    file MCFLemonSolver.h at line 339
 *  
 *  - V, which is the type of flows / deficits; typically, double can be used
 *    for maximum compatibility, but int (or even smaller) would yeld better
 *    performances;
 *  
 * - C, which is the type of ar costs; typically, double can be used for
 *    maximum compatibility, but int (or even smaller) would yeld better
 *    performances;
*/
template< typename GR, typename V, typename C>
class MCFLemonSolver<SMSppCapacityScaling, GR, V, C> : virtual public CDASolver, public MCFLemonSolverBase<SMSppCapacityScaling, GR, V, C>, Fields< SMSppCapacityScaling<GR, V, C> >
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

   

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
   kUnEval = 0     compute() has not been called yet

   kUnbounded = kUnEval + 1     the model is provably unbounded
 *  @{ */
  const int kErrorStatus = -1;
/*
    
    kUnEval = 0     compute() has not been called yet

    kUnbounded = kUnEval + 1     the model is provably unbounded

    kInfeasible                  the model is provably infeasible

    kBothInfeasible = kInfeasible + 1     both primal and dual infeasible

    kOK = 7         successful compute()
                    Any return value between kUnEval (excluded) and kOK
        (included) means that the object ran smoothly

    kStopTime = kOK + 1          stopped because of time limit

    kStopIter                    stopped because of iteration limit

    kError = 15     compute() stopped because of unrecoverable error
                    Any return value >= kError means that the object was
         forced to stop due to some error, e.g. of numerical nature

    kLowPrecision = kError + 1   a solution found but not provably optimal
    */



/** @} ---------------------------------------------------------------------*/
/*-------------- CONSTRUCTING AND DESTRUCTING MCFLemonSolver ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing MCFLemonSolver
 *  @{ */

 /// constructor: does nothing
 MCFLemonSolver( void ) : CDASolver(), MCFLemonSolverBase<SMSppCapacityScaling, GR, V, C>(), Fields<SMSppCapacityScaling<GR, V, C > >() {

 }
 
  /// destructor: doesnothing  
  ~MCFLemonSolver( void ) {   
  }
 
 

/*--------------------------------------------------------------------------*/
 /* intMaxIter = 0     maximum iterations for the next call to solve()

    intMaxSol          maximum number of different solutions to report

    intLogVerb         "verbosity" of the log

    intMaxDSol         maximum number of different dual solutions

    intLastParCDAS     first allowed parameter value for derived classes
    */

  enum LEMON_CS_int_par_type{
    intLastParLEMON_CS //< first allowed parameter value for derived classes
                           /**< convenience value for easily allow derived classes
                            * to further extend the set of types of return codes */
  };


/*--------------------------------------------------------------------------*/
 /* dblMaxTime = 0    maximum time for the next call to solve()

    dblRelAcc         relative accuracy for declaring a solution optimal

    dblAbsAcc          absolute accuracy for declaring a solution optimal

    dblUpCutOff        upper cutoff for stopping the algorithm

    dblLwCutOff        lower cutoff for stopping the algorithm

    dblRAccSol          maximum relative error in any reported solution

    dblAAccSol          maximum absolute error in any reported solution

    dblFAccSol          maximum constraint violation in any reported solution

    dblRAccDSol         maximum relative error in any dual solution

    dblAAccDSol         maximum absolute error in any dual solution

    dblFAccDSol         maximum absolute error in any dual solution

    dblLastParCDAS      first allowed parameter value for derived classes
    */

/*--------------------------------------------------------------------------*/

enum LEMON_CS_dbl_par_type{
        dblLastParLEMON_CS //< first allowed parameter value for derived classes
                           /**< convenience value for easily allow derived classes
                            * to further extend the set of types of return codes */
};



/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *
 * Parameter-wise, MCFSolver maps the parameters of [CDA]Solver
 *
 *  intMaxIter = 0    maximum iterations for the next call to solve()
 *  intMaxSol         maximum number of different solutions to report
 *  intLogVerb        "verbosity" of the log
 *  intMaxDSol        maximum number of different dual solutions
 *
 *  dblMaxTime = 0    maximum time for the next call to solve()
 *  dblRelAcc         relative accuracy for declaring a solution optimal
 *  dblAbsAcc         absolute accuracy for declaring a solution optimal
 *  dblUpCutOff       upper cutoff for stopping the algorithm
 *  dblLwCutOff       lower cutoff for stopping the algorithm
 *  dblRAccSol        maximum relative error in any reported solution
 *  dblAAccSol        maximum absolute error in any reported solution
 *  dblFAccSol        maximum constraint violation in any reported solution
 *  dblRAccDSol       maximum relative error in any dual solution
 *  dblAAccDSol       maximum absolute error in any dual solution
 *  dblFAccDSol       maximum absolute error in any dual solution
 *
 * into the parameter of MCFClass
 *
 * kMaxTime = 0       max time
 * kMaxIter           max number of iteration
 * kEpsFlw            tolerance for flows
 * kEpsDfct           tolerance for deficits
 * kEpsCst            tolerance for costs
 *
 * It then "extends" them, using
 *
 *  intLastParCDAS    first allowed parameter value for derived classes
 *  dblLastParCDAS    first allowed parameter value for derived classes
 *
 * In particular, one now has
 *
 * intLastParCDAS ==> kReopt             whether or not to reoptimize
 *
 * and any other parameter of specific :MCFClass following. This is done
 * via the two const static arrays Solver_2_MCFClass_int and
 * Solver_2_MCFClass_dbl, with a negative entry meaning "there is no such
 * parameter in MCFSolver".
 *
 *  @{ */



    /*--------------------------------------------------------------------------*/
    // set the ostream for the Solver log
    // not really, MCFClass objects are remarkably silent
    //
    // virtual void set_log( std::ostream *log_stream = nullptr ) override;

    /*--------------------------------------------------------------------------*/
  
    void set_par(idx_type par, int value){

        CDASolver::set_par(par,value);
        return;
        
    }
   
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    //TODO: change MCFC function to Algo function.
    /*void set_par(idx_type par, double value) override
    {
      if (Solver_2_MCFClass_dbl[par] >= 0)
      ;
      //  MCFC::SetPar(Solver_2_MCFClass_dbl[par], double(value));
    }*/
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/


    /** @} ---------------------------------------------------------------------*/
    /*--------------------- METHODS FOR SOLVING THE Block ----------------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Solving the MCF encoded by the current MCFBlock
     *  @{ */
    /// (try to) solve the MCF encoded in the MCFBlock 
        int compute( bool changedvars = true ) override {
        const static std::array<int, 6> LemonStatus_2_MCFstatus = {
                kErrorStatus, SMSpp_di_unipi_it::LEMON_sol_type::OPTIMAL,
                 kErrorStatus , SMSpp_di_unipi_it::LEMON_sol_type::INFEASIBLE,
                SMSpp_di_unipi_it::LEMON_sol_type::UNBOUNDED, kErrorStatus};
        
        const static std::array<int, 6> MCFstatus_2_sol_type = {
                kUnEval, Solver::kOK, kStopTime, kInfeasible, Solver::kUnbounded,
                Solver::kError};

        lock(); // first of all, acquire self-lock

        if (!f_Block)            // there is no [MCFBlock] to solve
                return (kBlockLocked); // return error

        bool owned = f_Block->is_owned_by(f_id); // check if already locked
        if ((!owned) && (!f_Block->read_lock())) // if not try to read_lock
                return (kBlockLocked);                 // return error on failure

        // while [read_]locked, process any outstanding Modification
        //TODO: ensure that modification are actually processed for MCFLemonSolver.
        //process_outstanding_Modification();

        if (!f_dmx_file.empty())
        { // if so required
                // output the current instance (after the changes) to a DMX file
                std::ofstream ProbFile(f_dmx_file, ios_base::out | ios_base::trunc);
                if (!ProbFile.is_open())
                throw(std::logic_error("cannot open DMX file " + f_dmx_file));

                //WriteMCF(ProbFile);
                writeDimacsMat(ProbFile, *MCFLemonSolver::dgp);
                ProbFile.close();
        }

        if (!owned)               // if the [MCF]Block was actually read_locked
                f_Block->read_unlock(); // read_unlock it

        auto start = chrono::system_clock::now();
        
      
        this->status = MCFLemonSolver::f_algo->run();

        
        auto end = chrono::system_clock::now();

        chrono::duration< double > elapsed = end - start;
        ticks = elapsed.count();

        
        
        unlock(); // release self-lock

        // now give out the result: note that the vector MCFstatus_2_sol_type[]
        // starts from 0 whereas the first value of MCFStatus is -1 (= kUnSolved),
        // hence the returned status has to be shifted by + 1
        return (MCFstatus_2_sol_type[this->get_status()]);
        }

    

    /** @} ---------------------------------------------------------------------*/
    /*---------------------- METHODS FOR READING RESULTS -----------------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Accessing the found solutions (if any)
     *  @{ */
    
    

    /// @brief getter for status field
    /// @param  
    /// @return value of this->status
    int get_status(void) const 
    {
      return (this->status);
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    ///Get elapsed time in run() method
    double get_elapsed_time(void) const override
    {
      return (this->ticks);
    }



    /*--------------------------------------------------------------------------*/
    //Return the lower bound solution(optimal) for the problem
    OFValue get_lb(void) override { return OFValue(MCFLemonSolver::f_algo->totalCost()); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    //Return the upper bound solution(optimal) for the problem
    OFValue get_ub(void) override { return OFValue(MCFLemonSolver::f_algo->totalCost()); }

    /*--------------------------------------------------------------------------*/
    //TODO: change MCFC function to Algo function. DONE
    bool has_var_solution(void) override
    {
     switch (this->get_status()) {
      case( SMSppCapacityScaling< GR , V , C >::ProblemType::OPTIMAL ):
      case( SMSppCapacityScaling< GR , V , C >::ProblemType::UNBOUNDED ):  
       return( true );
      default:
       return( false );
      }
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    //TODO: change MCFC function to Algo function.
    bool has_dual_solution(void) override
    {
     switch( this->get_status() ) {
      case( SMSppCapacityScaling< GR , V , C >::ProblemType::OPTIMAL ):
      case( SMSppCapacityScaling< GR , V , C >::ProblemType::INFEASIBLE ):
       return( true );
      default:
       return( false );
      }
    }

    /*--------------------------------------------------------------------------*/
    /*
     virtual bool is_var_feasible( void ) override { return( true ); }

     virtual bool is_dual_feasible( void ) override { return( true ); }
    */
    /*--------------------------------------------------------------------------*/
    /// write the "current" flow in the x ColVariable of the MCFBlock
    /** Write the "current" flow in the x ColVariable of the MCFBlock. To keep
     * the same format as MCFBlock::get_Solution() and
     * MCFBlock::map[forward/back]_Modification(), the Configuration *solc can
     * be used to "partly" save it. In particular, if solc != nullptr, it is
     * a SimpleConfiguration< int >, and solc->f_value == 2, then *nothing is
     * done*, since the Configuration is meant to say "only save/map the dual
     * solution". In all other cases, the flow solution is saved. */

    void get_var_solution(Configuration *solc = nullptr) override
    {/*
      if (!f_Block) // no [MCF]Block to write to
        return;     // cowardly and silently return

      auto tsolc = dynamic_cast<SimpleConfiguration<int> *>(solc);
      if (tsolc && (tsolc->f_value == 2))
        return;

      auto MCFB = static_cast<MCFBlock *>(f_Block);
      MCFBlock::Vec_FNumber X(MCFB->get_NArcs());
      this->MCFGetX(X.data());
      MCFB->set_x(X.begin());
   */}
  
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    /// write the "current" dual solution in the Constraint of the MCFBlock
    /** Write the "current" dual solution, i.e., node potentials and flow
     * reduced costs, in the dual variables of the Constraint (respectively,
     * the flow conservation constraints and bound ones) of the MCFBlock. To
     * keep the same format as MCFBlock::get_Solution() and
     * MCFBlock::map[forward/back]_Modification(), the Configuration *solc can
     * be used to "partly" save it. In particular, if solc != nullptr, it is
     * a SimpleConfiguration< int >, and solc->f_value == 1, then *nothing is
     * done*, since the Configuration is meant to say "only save/map the primal
     * solution". In all other cases, the flow solution is saved. */
    
    void get_dual_solution(Configuration *solc = nullptr) override
    {
    /*  if (!f_Block) // no [MCF]Block to write to
        return;     // cowardly and silently return

      auto tsolc = dynamic_cast<SimpleConfiguration<int> *>(solc);
      if (tsolc && (tsolc->f_value == 1))
        return;

      auto MCFB = static_cast<MCFBlock *>(f_Block);
      MCFBlock::Vec_CNumber Pi(MCFB->get_NNodes());
      this->MCFGetPi(Pi.data());
      MCFB->set_pi(Pi.begin());

      MCFBlock::Vec_FNumber RC(MCFB->get_NArcs());
      this->MCFGetRC(RC.data());
      MCFB->set_rc(RC.begin());
    */}
    
    /*--------------------------------------------------------------------------*/

   // bool new_var_solution(void) override { return (this->HaveNewX()); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   // bool new_dual_solution(void) override { return (this->HaveNewPi()); }

    /*--------------------------------------------------------------------------*/
    /*
     virtual void set_unbounded_threshold( const OFValue thr ) override { }
    */

    /*--------------------------------------------------------------------------*/

    bool has_var_direction(void) override { return (true); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    bool has_dual_direction(void) override { return (true); }

    /*--------------------------------------------------------------------------*/
    /// write the current direction in the x ColVariable of the MCFBlock
    /** Write the unbounded primal direction, i.e., augmenting cycle with
     * negative cost and unbounded capacity, in the x ColVariable of the
     * MCFBlock. To keep the same format as MCFBlock::get_Solution() and
     * MCFBlock::map[forward/back]_Modification(), the Configuration *solc can
     * be used to "partly" save it. In particular, if solc != nullptr, it is
     * a SimpleConfiguration< int >, and solc->f_value == 2, then *nothing is done*,
     * since the Configuration is meant to say "only save/map the dual
     * information". In all other cases, the direction (cycle) is saved.
     *
     * Or, rather, THIS SHOULD BE DONE, BUT THE METHOD IS NOT IMPLEMENTED yet. */

    void get_var_direction(Configuration *dirc = nullptr) override
    {
      auto tsolc = dynamic_cast<SimpleConfiguration<int> *>(dirc);
      if (tsolc && (tsolc->f_value == 2))
        return;

      throw(std::logic_error(
          "MCFSolver::get_var_direction() not implemented yet"));

      // TODO: implement using MCFC::MCFGetUnbCycl()
      // anyway, unsure if any current :MCFClass properly implemente the latter
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    /// write the current dual direction in the Constraint of the MCFBlock
    /** Write the current unbounded dual direction, i.e., a cut separating two
     * shores to that the residual demand in one is greater than the capacity
     * across them, in the Constraint of the Block, in particular in the dual
     * variables of the flow conservation ones. To keep the same format as
     * MCFBlock::get_Solution() and MCFBlock::map[forward/back]_Modification(),
     * the Configuration *solc can be used to "partly" save it. In particular, if
     * solc != nullptr, it is a SimpleConfiguration< int >, and solc->f_value == 1,
     * then *nothing is done*, since the Configuration is meant to say "only
     * save/map the primal information". In all other cases, the direction (cut)
     * is saved.
     *
     * Or, rather, THIS SHOULD BE DONE, BUT THE METHOD IS NOT IMPLEMENTED yet. */

    void get_dual_direction(Configuration *dirc = nullptr) override
    {
      auto tsolc = dynamic_cast<SimpleConfiguration<int> *>(dirc);
      if (tsolc && (tsolc->f_value == 1))
        return;

      throw(std::logic_error(
          "MCFSolver::get_dual_direction() not implemented yet"));

      // TODO: implement using MCFC::MCFGetUnfCut()
      // anyway, unsure if any current :MCFClass properly implemente the latter
    }

    /*--------------------------------------------------------------------------*/
    /*
     virtual bool new_var_direction( void ) override { return( false ); }

     virtual bool new_dual_direction( void ) override{ return( false ); }
    */
    /** @} ---------------------------------------------------------------------*/
    /*-------------- METHODS FOR READING THE DATA OF THE Solver ----------------*/
    /*--------------------------------------------------------------------------*/

    /*
     virtual bool is_dual_exact( void ) const override { return( true ); }
    */

    /*--------------------------------------------------------------------------*/
    /// "publicize" MCFClass::WriteMCF
    /** Make the method
     *
     *      void WriteMCF( ostream &oStrm , int frmt = 0 )
     *
     * of the base (private) MCFClass public, so that it can be freely used. */
    //TODO: change MCFC function to Algo function.
    //using MCFC::WriteMCF;

    /*--------------------------------------------------------------------------*/
    /*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Handling the parameters of the MCFLemonSolver
     *
     * Each MCFLemonSolver< Algo > may have its own extra int / double parameters. If
     * this is the case, it will have to specialize the following methods to
     * handle them. The general definition just handles the case of the
     *
     * intLastParCDAS ==> kReopt             whether or not to reoptimize
     *
     *  @{ */

     ///@return number of int algorithmic parameters
    [[nodiscard]] idx_type get_num_int_par(void) const override
    {
        return(intLastParLEMON_CS);
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @return number of dbl algorithmic parameters
    [[nodiscard]] idx_type get_num_dbl_par(void) const override
    {
        return(dblLastParLEMON_CS);
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @return number of str algorithimc parameters 
    [[nodiscard]] idx_type get_num_str_par(void) const override
    {
      return (MCFLemonSolver::strLastParLEMON);
    }
    
    /*--------------------------------------------------------------------------*/
    
    /// @brief used for obtain default value of an int parameter
    /// @param par 
    /// @return default value of parameter par, if exists
    [[nodiscard]] int get_dflt_int_par(idx_type par) const override
    {
      if(par > intLastParLEMON_CS){
        throw std::invalid_argument(std::to_string(par));
      }

      switch(par){
        default: return (CDASolver::get_dflt_int_par(par));
      }
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for obtain default value of an dbl parameter
    /// @param par 
    /// @return default value of parameter par, if exists
    [[nodiscard]] double get_dflt_dbl_par(idx_type par) const override
    {
       if(par > dblLastParLEMON_CS){
          throw std::invalid_argument(std::to_string(par));
       }    
       
       switch(par){
        default: return (CDASolver::get_dflt_dbl_par(par));
       }
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for obtain default value of an string parameter
    /// @param par 
    /// @return default value of parameter par, if exists
    [[nodiscard]] const std::string &get_dflt_str_par(idx_type par)
        const override
    {

      if(par > MCFLemonSolverBase<SMSppCapacityScaling, GR, V, C>::strLastParLEMON){
        throw std::invalid_argument("Invalid str parameter: out_of_range " + std::to_string(par));
      }

      static const std::string _empty;
      if (par == MCFLemonSolverBase<SMSppCapacityScaling, GR, V, C>::strLastParLEMON)
        return (_empty);

      return (CDASolver::get_dflt_str_par(par));
    }
    
    /*--------------------------------------------------------------------------*/
    
    /// @brief used for get the value of int parameters
    /// @param par 
    /// @return value of parameter indexed by par
    [[nodiscard]] int get_int_par(idx_type par) const override
    {
      //No int parameters
      return (get_dflt_int_par(par));
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for get the value of dbl parameters
    /// @param par 
    /// @return value of parameter indexed by par
    [[nodiscard]] double get_dbl_par(idx_type par) const override
    {
      //No dbl parameters
      return (get_dflt_dbl_par(par));
    }
        
        
    /*--------------------------------------------------------------------------*/
    
    /// @brief used for get the value of string parameters
    /// @param par 
    /// @return value of parameter indexed by par
    [[nodiscard]] const std::string & get_str_par( idx_type par ) const override {
        
        if( par == MCFLemonSolverBase<SMSppCapacityScaling, GR, V, C>::strDMXFile )
        return( this->f_dmx_file );

        return( get_dflt_str_par( par ) );
        }
   

    /*--------------------------------------------------------------------------*/

    /// @brief used for convert string to int parameter's index
    /// @param name 
    /// @return an index that denotes parameter name
    [[nodiscard]] idx_type int_par_str2idx(const std::string &name)
        const override
    {
      return (CDASolver::int_par_str2idx(name));
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for convert string to dbl parameter's index
    /// @param name 
    /// @return an index that denotes parameter name
    [[nodiscard]] idx_type dbl_par_str2idx(const std::string &name)
        const override
    {
      return (CDASolver::dbl_par_str2idx(name));
    }
      
    /*--------------------------------------------------------------------------*/
    
    /// @brief used for convert int index idx to phrasal rapresentation
    /// @param idx 
    /// @return a string that denotes index idx parameter
    [[nodiscard]] const std::string &int_par_idx2str(idx_type idx)
        const override
    {
      if(idx > intLastParLEMON_CS){
        throw std::invalid_argument(std::to_string(idx));
      }

      switch(idx){
        default: break;
      }

      return (CDASolver::int_par_idx2str(idx));
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for convert dbl index idx to phrasal rapresentation
    /// @param idx 
    /// @return a string that denotes index idx parameter
    [[nodiscard]] const std::string &dbl_par_idx2str(idx_type idx)
        const override
    {
      if(idx > dblLastParLEMON_CS){
        throw std::invalid_argument(std::to_string(idx));
      }

      switch(idx){
        default: break;
      }

      return (CDASolver::dbl_par_idx2str(idx));
    }

    
    /** @} ---------------------------------------------------------------------*/
    /*------------ METHODS FOR HANDLING THE State OF THE MCFLemonSolver -------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Handling the State of the MCFSolver
     *  @{ */

    //[[nodiscard]] State *get_State(void) const override;

    /*--------------------------------------------------------------------------*/

    //void put_State(const State &state) override;

    /*--------------------------------------------------------------------------*/

    //void put_State(State &&state) override;

    /** @} ---------------------------------------------------------------------*/
    /*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Changing the data of the model
     *  @{ */

    /** The only reason why MCFSolver::add_Modification() needs be defined is to
     * properly react to NBModification. Indeed, the correct reaction is to
     * *immediately* reload the MCFBlock, besides clearing the list of
     * Modification as Solver::add_Modification() already does. The issue is
     * that if arcs/nodes are added/deleted after the NBModification is issued
     * but before it is processed, then the number of nodes/arcs at the moment
     * in which the NBModification is processed is different from that at the
     * moment in which is issued, which may break the "naming convention"
     * (because the name of, say, a newly created arc depends on the current
     * state and/or number of the arcs).
     *
     * Important note: THIS VERSION ONLY WORKS PROPERLY IF THE MCFBlock IS
     * "FRESHLY MINTED", I.E., THERE ARE NO CLOSED OR DELETED ARCS.
     *
     * This should ordinarily always happen, as whenever the MCFBlock is changed
     * the NBModification is immediately issued. The problem may come if the
     * MCFBlock is a R3Block of another MCFBlock which is loaded and then
     * further modified, and the NBModification to this MCFBlock is generated by
     * a map_forward_Modification() of the NBModification to the original
     * MCFBlock: then, this MCFBlock may be copied from a MCFBlock that has
     * closed or deleted arcs and this method would not work. */

    /*void add_Modification(sp_Mod &mod) override
    {
      if (std::dynamic_pointer_cast<const NBModification>(mod))
      {
        // this is the "nuclear option": the MCFBlock has been re-loaded, so
        // the MCFClass solver also has to (immediately)
        //TODO: change MCFC function to Algo function.
        auto MCFB = static_cast<MCFBlock *>(f_Block);
        MCFC::LoadNet(MCFB->get_MaxNNodes(), MCFB->get_MaxNArcs(),
                      MCFB->get_NNodes(), MCFB->get_NArcs(),
                      MCFB->get_U().empty() ? nullptr : MCFB->get_U().data(),
                      MCFB->get_C().empty() ? nullptr : MCFB->get_C().data(),
                      MCFB->get_B().empty() ? nullptr : MCFB->get_B().data(),
                      MCFB->get_SN().data(), MCFB->get_EN().data());
        // TODO: PreProcess() changes the internal data of the MCFSolver using
        //       information about how the data of the MCF is *now*. If the
        //       data changes, some of the deductions (say, reducing the capacity
        //       of some arcs) may no longer be correct and they should be undone,
        //       but there isn't any proper way to handle this. Thus, PreProcess()
        //       has to be disabled for now; maybe later on someone will take
        //       care to make this work (or maybe not).
        // MCFC::PreProcess();
        // besides, any outstanding modification makes no sense any longer
        mod_clear();
      }
      else
        push_back(mod);
    }*/

    /*--------------------------------------------------------------------------*/
    /*-------------------------------- FRIENDS ---------------------------------*/
    /*--------------------------------------------------------------------------*/

    friend class MCFLemonState; // make MCFSolverState friend

    /** @} ---------------------------------------------------------------------*/
    /*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
    /*--------------------------------------------------------------------------*/

  protected:
    /*--------------------------------------------------------------------------*/
    /*-------------------------- PROTECTED METHODS -----------------------------*/
    /*--------------------------------------------------------------------------*/

    void process_outstanding_Modification(void);

    void guts_of_poM(c_p_Mod mod);

    /*--------------------------------------------------------------------------*/
    /*---------------------------- PROTECTED FIELDS  ---------------------------*/
    /*--------------------------------------------------------------------------*/

    const static std::vector<int> Solver_2_MCFClass_int;
    // the (static const) map between Solver int parameters and MCFClass ones

    const static std::vector<int> Solver_2_MCFClass_dbl;
    // the (static const) map between Solver int parameters and MCFClass ones

    std::string f_dmx_file; 
    // string for DMX file output


    /*--------------------------------------------------------------------------*/
    /*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
    /*--------------------------------------------------------------------------*/

  private:
    /*--------------------------------------------------------------------------*/
    /*-------------------------- PRIVATE METHODS -------------------------------*/
    /*--------------------------------------------------------------------------*/

    SMSpp_insert_in_factory_h;

     /*--------------------------------------------------------------------------*/
    /*-------------------------- PRIVATE FIELDS -------------------------------*/
    /*--------------------------------------------------------------------------*/
    int status = SMSpp_di_unipi_it::LEMON_sol_type::UNSOLVED;  //Variable used in compute function for getting status
    typename SMSppCapacityScaling< GR , V , C >::ProblemType status_2_pType;
    //Status of compute() method
    
    double ticks;  //Elaped time in ticks for compute() method
    /*--------------------------------------------------------------------------*/

  }; // end( class MCFLemonSolver<SMSppCapacityScaling, GR, V, C>  Specialization)


/*--------------------------------------------------------------------------*/
/*-------------------------------- SPECIALIZED CLASSES ---------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------------- COSTSCALING -----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/


/** Specialized MCFLemonSolver<CostScaling, GR, V, C> that contains specialized
 *  compute() method, enums for indexing algorithimc parameters and function
 *  set/get_*_par for manage them.
 *  This class derives from MCFLemonSolverBase for using set_Block() methods and
 *  for inherit f_algo and dgp fields.
 * 
 *  Template parameters are:
 * 
 *   - CostScaling implements a cost scaling algorithm that performs push/augment and
 *     relabel operations for finding a minimum cost flow. It is a highly efficient primal-dual
 *     solution method, which can be viewed as the generalization of the preflow push-relabel
 *     algorithm for the maximum flow problem. It is a polynomial algorithm.
 * 
 *  - GR that represents the directed graph, the possibilities are described in the
 *    file MCFLemonSolver.h at line 339
 *  
 *  - V, which is the type of flows / deficits; typically, double can be used
 *    for maximum compatibility, but int (or even smaller) would yeld better
 *    performances;
 *  
 * - C, which is the type of ar costs; typically, double can be used for
 *    maximum compatibility, but int (or even smaller) would yeld better
 *    performances;
*/

template< typename GR, typename V, typename C>
class MCFLemonSolver<SMSppCostScaling, GR, V, C> : virtual public CDASolver, public MCFLemonSolverBase<SMSppCostScaling, GR, V, C>, Fields< SMSppCostScaling<GR, V, C> >
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

   

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
   kUnEval = 0     compute() has not been called yet

   kUnbounded = kUnEval + 1     the model is provably unbounded
 *  @{ */
  const int kErrorStatus = -1;
/*
    
    kUnEval = 0     compute() has not been called yet

    kUnbounded = kUnEval + 1     the model is provably unbounded

    kInfeasible                  the model is provably infeasible

    kBothInfeasible = kInfeasible + 1     both primal and dual infeasible

    kOK = 7         successful compute()
                    Any return value between kUnEval (excluded) and kOK
        (included) means that the object ran smoothly

    kStopTime = kOK + 1          stopped because of time limit

    kStopIter                    stopped because of iteration limit

    kError = 15     compute() stopped because of unrecoverable error
                    Any return value >= kError means that the object was
         forced to stop due to some error, e.g. of numerical nature

    kLowPrecision = kError + 1   a solution found but not provably optimal
    */



/** @} ---------------------------------------------------------------------*/
/*-------------- CONSTRUCTING AND DESTRUCTING MCFLemonSolver ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing MCFLemonSolver
 *  @{ */

 /// constructor: Build f_method_type
 /** Void constructor: Initialize f_method_type to default algorithmic parameter */

 MCFLemonSolver( void ) : CDASolver(), Fields<SMSppCostScaling<GR, V, C > >() {
  f_method_type = SMSppCostScaling<GR, V, C>::Method::PARTIAL_AUGMENT;


 /* static_assert( std::is_same< Algo< GR , V , C > ,
                               SMSppCapacityScaling< GR , V , C > >::value ||
		 std::is_same< Algo< GR , V , C > ,
                               SMSppCostScaling< GR , V , C >::value     ||
		 std::is_same< Algo< GR , V , C > ,
                               CycleCanceling< GR , V , C > >::value       ||
		 std::is_same< Algo< GR , V , C > ,
		               NetworkSimplex< GR , V , C > >::value ,
		 "Algo must be one of the LEMON algorithms");
  */
  }
  
  /// @brief delete f_method algorithmic parameter pointer
  /// @param  
  ~MCFLemonSolver( void ) {
    delete Fields<SMSppCostScaling<GR, V, C>>::f_method;    
  }
 
 

/*--------------------------------------------------------------------------*/
 /* intMaxIter = 0     maximum iterations for the next call to solve()

    intMaxSol          maximum number of different solutions to report

    intLogVerb         "verbosity" of the log

    intMaxDSol         maximum number of different dual solutions

    intLastParCDAS     first allowed parameter value for derived classes
    */

   enum LEMON_CS_int_par_type{
        kMethod = intLastParCDAS,
        intLastParLEMON_CS //< first allowed parameter value for derived classes
                           /**< convenience value for easily allow derived classes
                            * to further extend the set of types of return codes */
   };

/*--------------------------------------------------------------------------*/
 /* dblMaxTime = 0    maximum time for the next call to solve()

    dblRelAcc         relative accuracy for declaring a solution optimal

    dblAbsAcc          absolute accuracy for declaring a solution optimal

    dblUpCutOff        upper cutoff for stopping the algorithm

    dblLwCutOff        lower cutoff for stopping the algorithm

    dblRAccSol          maximum relative error in any reported solution

    dblAAccSol          maximum absolute error in any reported solution

    dblFAccSol          maximum constraint violation in any reported solution

    dblRAccDSol         maximum relative error in any dual solution

    dblAAccDSol         maximum absolute error in any dual solution

    dblFAccDSol         maximum absolute error in any dual solution

    dblLastParCDAS      first allowed parameter value for derived classes
    */

/*--------------------------------------------------------------------------*/
 enum LEMON_CS_dbl_par_type{
        dblLastParLEMON_CS //< first allowed parameter value for derived classes
                           /**< convenience value for easily allow derived classes
                            * to further extend the set of types of return codes */
 };

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *
 * Parameter-wise, MCFSolver maps the parameters of [CDA]Solver
 *
 *  intMaxIter = 0    maximum iterations for the next call to solve()
 *  intMaxSol         maximum number of different solutions to report
 *  intLogVerb        "verbosity" of the log
 *  intMaxDSol        maximum number of different dual solutions
 *
 *  dblMaxTime = 0    maximum time for the next call to solve()
 *  dblRelAcc         relative accuracy for declaring a solution optimal
 *  dblAbsAcc         absolute accuracy for declaring a solution optimal
 *  dblUpCutOff       upper cutoff for stopping the algorithm
 *  dblLwCutOff       lower cutoff for stopping the algorithm
 *  dblRAccSol        maximum relative error in any reported solution
 *  dblAAccSol        maximum absolute error in any reported solution
 *  dblFAccSol        maximum constraint violation in any reported solution
 *  dblRAccDSol       maximum relative error in any dual solution
 *  dblAAccDSol       maximum absolute error in any dual solution
 *  dblFAccDSol       maximum absolute error in any dual solution
 *
 * into the parameter of MCFClass
 *
 * kMaxTime = 0       max time
 * kMaxIter           max number of iteration
 * kEpsFlw            tolerance for flows
 * kEpsDfct           tolerance for deficits
 * kEpsCst            tolerance for costs
 *
 * It then "extends" them, using
 *
 *  intLastParCDAS    first allowed parameter value for derived classes
 *  dblLastParCDAS    first allowed parameter value for derived classes
 *
 * In particular, one now has
 *
 * intLastParCDAS ==> kReopt             whether or not to reoptimize
 *
 * and any other parameter of specific :MCFClass following. This is done
 * via the two const static arrays Solver_2_MCFClass_int and
 * Solver_2_MCFClass_dbl, with a negative entry meaning "there is no such
 * parameter in MCFSolver".
 *
 *  @{ */



    /*--------------------------------------------------------------------------*/
    // set the ostream for the Solver log
    // not really, MCFClass objects are remarkably silent
    //
    // virtual void set_log( std::ostream *log_stream = nullptr ) override;

    /*--------------------------------------------------------------------------*/
    
    /// @brief set the parameter par with value
    /// @param par 
    /// @param value 
    void set_par(idx_type par, int value) 
    {
        if(par != kMethod){
                CDASolver::set_par(par, value);
                return;
        }
        
        if((value < 0) || (value > 4)){
                throw(std::invalid_argument(std::to_string(value)));
        }

        f_method_type = value;
        free(Fields<SMSppCostScaling<GR, V, C>>::f_method);
        Fields<SMSppCostScaling<GR, V, C>>::f_method = NULL;
        return;

      
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    //TODO: change MCFC function to Algo function.
    /*void set_par(idx_type par, double value) override
    {
      if (Solver_2_MCFClass_dbl[par] >= 0)
      ;
      //  MCFC::SetPar(Solver_2_MCFClass_dbl[par], double(value));
    }*/
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/


    /** @} ---------------------------------------------------------------------*/
    /*--------------------- METHODS FOR SOLVING THE Block ----------------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Solving the MCF encoded by the current MCFBlock
     *  @{ */
    /// (try to) solve the MCF encoded in the MCFBlock 
        int compute( bool changedvars = true ) override {
        const static std::array<int, 6> LemonStatus_2_MCFstatus = {
                kErrorStatus, SMSpp_di_unipi_it::LEMON_sol_type::OPTIMAL,
                 kErrorStatus , SMSpp_di_unipi_it::LEMON_sol_type::INFEASIBLE,
                SMSpp_di_unipi_it::LEMON_sol_type::UNBOUNDED, kErrorStatus};
        
        const static std::array<int, 6> MCFstatus_2_sol_type = {
                kUnEval, Solver::kOK, kStopTime, kInfeasible, Solver::kUnbounded,
                Solver::kError};

        lock(); // first of all, acquire self-lock

        if (!f_Block)            // there is no [MCFBlock] to solve
                return (kBlockLocked); // return error

        bool owned = f_Block->is_owned_by(f_id); // check if already locked
        if ((!owned) && (!f_Block->read_lock())) // if not try to read_lock
                return (kBlockLocked);                 // return error on failure

        // while [read_]locked, process any outstanding Modification
        //TODO: ensure that modification are actually processed for MCFLemonSolver.
        //process_outstanding_Modification();

        if (!f_dmx_file.empty())
        { // if so required
                // output the current instance (after the changes) to a DMX file
                std::ofstream ProbFile(f_dmx_file, ios_base::out | ios_base::trunc);
                if (!ProbFile.is_open())
                throw(std::logic_error("cannot open DMX file " + f_dmx_file));

                //WriteMCF(ProbFile);
                writeDimacsMat(ProbFile, *MCFLemonSolver::dgp);
                ProbFile.close();
        }

        if (!owned)               // if the [MCF]Block was actually read_locked
                f_Block->read_unlock(); // read_unlock it

        auto start = chrono::system_clock::now();
        
        if(*Fields<SMSppCostScaling< GR, V, C > >::f_method != NULL){
        this->status = MCFLemonSolver::f_algo->run(*Fields<SMSppCostScaling< GR, V, C > >::f_method);
        }else{
                //Build f_method and execute run() method
                Fields<SMSppCostScaling<GR, V, C> >::f_method = new typename SMSppCostScaling<GR, V, C>::Method(f_method_type);
                this->status = MCFLemonSolver::f_algo->run(*Fields<SMSppCostScaling< GR, V, C > >::f_method);
                
        }
        
        auto end = chrono::system_clock::now();

        chrono::duration< double > elapsed = end - start;
        ticks = elapsed.count();

        
        
        unlock(); // release self-lock

        // now give out the result: note that the vector MCFstatus_2_sol_type[]
        // starts from 0 whereas the first value of MCFStatus is -1 (= kUnSolved),
        // hence the returned status has to be shifted by + 1
        return (MCFstatus_2_sol_type[this->get_status()]);
        }

    

    /** @} ---------------------------------------------------------------------*/
    /*---------------------- METHODS FOR READING RESULTS -----------------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Accessing the found solutions (if any)
     *  @{ */
    
    

    /// @brief getter of the status
    /// @param  
    /// @return this->status
    int get_status(void) const 
    {
      return (this->status);
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief getter of ticks (elapsed time of run())
    ///@return this->ticks
    double get_elapsed_time(void) const override
    {
      return (this->ticks);
    }



    /*--------------------------------------------------------------------------*/
    //Return the lower bound solution(optimal) for the problem
    OFValue get_lb(void) override { return OFValue(MCFLemonSolver::f_algo->totalCost()); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    //Return the upper bound solution(optimal) for the problem
    OFValue get_ub(void) override { return OFValue(MCFLemonSolver::f_algo->totalCost()); }

    /*--------------------------------------------------------------------------*/
    //TODO: change MCFC function to Algo function. DONE
    bool has_var_solution(void) override
    {
     switch (this->get_status()) {
      case( SMSppCostScaling< GR , V , C >::ProblemType::OPTIMAL ):
      case( SMSppCostScaling< GR , V , C >::ProblemType::UNBOUNDED ):  
       return( true );
      default:
       return( false );
      }
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    //TODO: change MCFC function to Algo function.
    bool has_dual_solution(void) override
    {
     switch( this->get_status() ) {
      case( SMSppCostScaling< GR , V , C >::ProblemType::OPTIMAL ):
      case( SMSppCostScaling< GR , V , C >::ProblemType::INFEASIBLE ):
       return( true );
      default:
       return( false );
      }
    }

    /*--------------------------------------------------------------------------*/
    /*
     virtual bool is_var_feasible( void ) override { return( true ); }

     virtual bool is_dual_feasible( void ) override { return( true ); }
    */
    /*--------------------------------------------------------------------------*/
    /// write the "current" flow in the x ColVariable of the MCFBlock
    /** Write the "current" flow in the x ColVariable of the MCFBlock. To keep
     * the same format as MCFBlock::get_Solution() and
     * MCFBlock::map[forward/back]_Modification(), the Configuration *solc can
     * be used to "partly" save it. In particular, if solc != nullptr, it is
     * a SimpleConfiguration< int >, and solc->f_value == 2, then *nothing is
     * done*, since the Configuration is meant to say "only save/map the dual
     * solution". In all other cases, the flow solution is saved. */

    void get_var_solution(Configuration *solc = nullptr) override
    {/*
      if (!f_Block) // no [MCF]Block to write to
        return;     // cowardly and silently return

      auto tsolc = dynamic_cast<SimpleConfiguration<int> *>(solc);
      if (tsolc && (tsolc->f_value == 2))
        return;

      auto MCFB = static_cast<MCFBlock *>(f_Block);
      MCFBlock::Vec_FNumber X(MCFB->get_NArcs());
      this->MCFGetX(X.data());
      MCFB->set_x(X.begin());
   */}
  
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    /// write the "current" dual solution in the Constraint of the MCFBlock
    /** Write the "current" dual solution, i.e., node potentials and flow
     * reduced costs, in the dual variables of the Constraint (respectively,
     * the flow conservation constraints and bound ones) of the MCFBlock. To
     * keep the same format as MCFBlock::get_Solution() and
     * MCFBlock::map[forward/back]_Modification(), the Configuration *solc can
     * be used to "partly" save it. In particular, if solc != nullptr, it is
     * a SimpleConfiguration< int >, and solc->f_value == 1, then *nothing is
     * done*, since the Configuration is meant to say "only save/map the primal
     * solution". In all other cases, the flow solution is saved. */
    
    void get_dual_solution(Configuration *solc = nullptr) override
    {
    /*  if (!f_Block) // no [MCF]Block to write to
        return;     // cowardly and silently return

      auto tsolc = dynamic_cast<SimpleConfiguration<int> *>(solc);
      if (tsolc && (tsolc->f_value == 1))
        return;

      auto MCFB = static_cast<MCFBlock *>(f_Block);
      MCFBlock::Vec_CNumber Pi(MCFB->get_NNodes());
      this->MCFGetPi(Pi.data());
      MCFB->set_pi(Pi.begin());

      MCFBlock::Vec_FNumber RC(MCFB->get_NArcs());
      this->MCFGetRC(RC.data());
      MCFB->set_rc(RC.begin());
    */}
    
    /*--------------------------------------------------------------------------*/

   // bool new_var_solution(void) override { return (this->HaveNewX()); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   // bool new_dual_solution(void) override { return (this->HaveNewPi()); }

    /*--------------------------------------------------------------------------*/
    /*
     virtual void set_unbounded_threshold( const OFValue thr ) override { }
    */

    /*--------------------------------------------------------------------------*/

    bool has_var_direction(void) override { return (true); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    bool has_dual_direction(void) override { return (true); }

    /*--------------------------------------------------------------------------*/
    /// write the current direction in the x ColVariable of the MCFBlock
    /** Write the unbounded primal direction, i.e., augmenting cycle with
     * negative cost and unbounded capacity, in the x ColVariable of the
     * MCFBlock. To keep the same format as MCFBlock::get_Solution() and
     * MCFBlock::map[forward/back]_Modification(), the Configuration *solc can
     * be used to "partly" save it. In particular, if solc != nullptr, it is
     * a SimpleConfiguration< int >, and solc->f_value == 2, then *nothing is done*,
     * since the Configuration is meant to say "only save/map the dual
     * information". In all other cases, the direction (cycle) is saved.
     *
     * Or, rather, THIS SHOULD BE DONE, BUT THE METHOD IS NOT IMPLEMENTED yet. */

    void get_var_direction(Configuration *dirc = nullptr) override
    {
      auto tsolc = dynamic_cast<SimpleConfiguration<int> *>(dirc);
      if (tsolc && (tsolc->f_value == 2))
        return;

      throw(std::logic_error(
          "MCFSolver::get_var_direction() not implemented yet"));

      // TODO: implement using MCFC::MCFGetUnbCycl()
      // anyway, unsure if any current :MCFClass properly implemente the latter
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    /// write the current dual direction in the Constraint of the MCFBlock
    /** Write the current unbounded dual direction, i.e., a cut separating two
     * shores to that the residual demand in one is greater than the capacity
     * across them, in the Constraint of the Block, in particular in the dual
     * variables of the flow conservation ones. To keep the same format as
     * MCFBlock::get_Solution() and MCFBlock::map[forward/back]_Modification(),
     * the Configuration *solc can be used to "partly" save it. In particular, if
     * solc != nullptr, it is a SimpleConfiguration< int >, and solc->f_value == 1,
     * then *nothing is done*, since the Configuration is meant to say "only
     * save/map the primal information". In all other cases, the direction (cut)
     * is saved.
     *
     * Or, rather, THIS SHOULD BE DONE, BUT THE METHOD IS NOT IMPLEMENTED yet. */

    void get_dual_direction(Configuration *dirc = nullptr) override
    {
      auto tsolc = dynamic_cast<SimpleConfiguration<int> *>(dirc);
      if (tsolc && (tsolc->f_value == 1))
        return;

      throw(std::logic_error(
          "MCFSolver::get_dual_direction() not implemented yet"));

      // TODO: implement using MCFC::MCFGetUnfCut()
      // anyway, unsure if any current :MCFClass properly implemente the latter
    }

    /*--------------------------------------------------------------------------*/
    /*
     virtual bool new_var_direction( void ) override { return( false ); }

     virtual bool new_dual_direction( void ) override{ return( false ); }
    */
    /** @} ---------------------------------------------------------------------*/
    /*-------------- METHODS FOR READING THE DATA OF THE Solver ----------------*/
    /*--------------------------------------------------------------------------*/

    /*
     virtual bool is_dual_exact( void ) const override { return( true ); }
    */

    /*--------------------------------------------------------------------------*/
    /// "publicize" MCFClass::WriteMCF
    /** Make the method
     *
     *      void WriteMCF( ostream &oStrm , int frmt = 0 )
     *
     * of the base (private) MCFClass public, so that it can be freely used. */
    //TODO: change MCFC function to Algo function.
    //using MCFC::WriteMCF;

    /*--------------------------------------------------------------------------*/
    /*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Handling the parameters of the MCFLemonSolver
     *
     * Each MCFLemonSolver< Algo > may have its own extra int / double parameters. If
     * this is the case, it will have to specialize the following methods to
     * handle them. The general definition just handles the case of the
     *
     * intLastParCDAS ==> kReopt             whether or not to reoptimize
     *
     *  @{ */

    ///@return number of int algorithmic parameters
    [[nodiscard]] idx_type get_num_int_par(void) const override
    { 
       return (intLastParLEMON_CS);
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    ///@return number of dbl algorithmic parameters
    [[nodiscard]] idx_type get_num_dbl_par(void) const override
    {
        return (dblLastParLEMON_CS);
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    ///@return number of str algorithmic parameters
    [[nodiscard]] idx_type get_num_str_par(void) const override
    {
      return (MCFLemonSolver::strLastParLEMON);
    }
    
    /*--------------------------------------------------------------------------*/
    
    /// @brief used for obtain default value of an int parameter
    /// @param par 
    /// @return default value of parameter par, if exists
    [[nodiscard]] int get_dflt_int_par(idx_type par) const override
    {

      if(par > intLastParLEMON_CS){
        throw std::invalid_argument(std::to_string(par));
      }

      switch(par){
        case kMethod: return SMSppCostScaling<GR, V, C>::Method::PARTIAL_AUGMENT;
        default: return (CDASolver::get_dflt_int_par(par));
      }
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for obtain default value of an dbl parameter
    /// @param par 
    /// @return default value of parameter par, if exists
    [[nodiscard]] double get_dflt_dbl_par(idx_type par) const override
    {
      if(par > intLastParLEMON_CS){
        throw std::invalid_argument(std::to_string(par));
      }
      
      switch(par){
        default: return(CDASolver::get_dflt_dbl_par(par));
      }
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for obtain default value of an string parameter
    /// @param par 
    /// @return default value of parameter par, if exists
    [[nodiscard]] const std::string &get_dflt_str_par(idx_type par)
        const override
    {

      if(par > MCFLemonSolverBase<SMSppCostScaling, GR, V, C>::strLastParLEMON){
        throw std::invalid_argument("Invalid str parameter: out_of_range " + std::to_string(par));
      }

      static const std::string _empty;
      if (par == MCFLemonSolverBase<SMSppCostScaling, GR, V, C>::strLastParLEMON)
        return (_empty);

      return (CDASolver::get_dflt_str_par(par));
    }
    
    /*--------------------------------------------------------------------------*/
    
    /// @brief used for get the value of int parameters
    /// @param par 
    /// @return value of parameter indexed by par
    [[nodiscard]] int get_int_par(idx_type par) const override
    {

      if(par == kMethod){
        return f_method_type;
      }
      return (get_dflt_int_par(par));
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for get the value of dbl parameters
    /// @param par 
    /// @return value of parameter indexed by par
    [[nodiscard]] double get_dbl_par(idx_type par) const override
    {
      return (get_dflt_dbl_par(par));
    }
        
    /*--------------------------------------------------------------------------*/
    
    /// @brief used for get the value of string parameters
    /// @param par 
    /// @return value of parameter indexed by par
    [[nodiscard]] const std::string & get_str_par( idx_type par ) const override {
        
        if( par == MCFLemonSolverBase<SMSppCostScaling, GR, V, C>::strDMXFile )
        return( this->f_dmx_file );

        return( get_dflt_str_par( par ) );
    }  
    
    
    /*--------------------------------------------------------------------------*/

    /// @brief used for convert string to int parameter's index
    /// @param name 
    /// @return an index that denotes parameter name
    [[nodiscard]] idx_type int_par_str2idx(const std::string &name)
        const override
    {
      if (name == "kMethod")
        return (kMethod);

      return (CDASolver::int_par_str2idx(name));
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for convert string to dbl parameter's index
    /// @param name 
    /// @return an index that denotes parameter name
    [[nodiscard]] idx_type dbl_par_str2idx(const std::string &name)
        const override
    {
      return ( name == "kMethod" ? kMethod : CDASolver::dbl_par_str2idx(name));
    }
      
    /*--------------------------------------------------------------------------*/
    
    /// @brief used for convert int index idx to phrasal rapresentation
    /// @param idx 
    /// @return a string that denotes index idx parameter
    [[nodiscard]] const std::string &int_par_idx2str(idx_type idx)
        const override
    {
      if(idx > intLastParLEMON_CS){
        throw std::invalid_argument(std::to_string(idx));
      }
      
      static const std::string par = "kMethod";
      switch(idx){
        case kMethod: return (par);
        default: break;
      }

      return (CDASolver::int_par_idx2str(idx));
    }
    
    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    
    /// @brief used for convert dbl index idx to phrasal rapresentation
    /// @param idx 
    /// @return a string that denotes index idx parameter
    [[nodiscard]] const std::string &dbl_par_idx2str(idx_type idx)
        const override
    {
      if(idx > dblLastParLEMON_CS){
        throw std::invalid_argument(std::to_string(idx));
      }

      switch(idx){
        default: break;
      }

      return (CDASolver::dbl_par_idx2str(idx));
    }
    
    /** @} ---------------------------------------------------------------------*/
    /*------------ METHODS FOR HANDLING THE State OF THE MCFLemonSolver -------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Handling the State of the MCFSolver
     *  @{ */

    //[[nodiscard]] State *get_State(void) const override;

    /*--------------------------------------------------------------------------*/

    //void put_State(const State &state) override;

    /*--------------------------------------------------------------------------*/

    //void put_State(State &&state) override;

    /** @} ---------------------------------------------------------------------*/
    /*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Changing the data of the model
     *  @{ */

    /** The only reason why MCFSolver::add_Modification() needs be defined is to
     * properly react to NBModification. Indeed, the correct reaction is to
     * *immediately* reload the MCFBlock, besides clearing the list of
     * Modification as Solver::add_Modification() already does. The issue is
     * that if arcs/nodes are added/deleted after the NBModification is issued
     * but before it is processed, then the number of nodes/arcs at the moment
     * in which the NBModification is processed is different from that at the
     * moment in which is issued, which may break the "naming convention"
     * (because the name of, say, a newly created arc depends on the current
     * state and/or number of the arcs).
     *
     * Important note: THIS VERSION ONLY WORKS PROPERLY IF THE MCFBlock IS
     * "FRESHLY MINTED", I.E., THERE ARE NO CLOSED OR DELETED ARCS.
     *
     * This should ordinarily always happen, as whenever the MCFBlock is changed
     * the NBModification is immediately issued. The problem may come if the
     * MCFBlock is a R3Block of another MCFBlock which is loaded and then
     * further modified, and the NBModification to this MCFBlock is generated by
     * a map_forward_Modification() of the NBModification to the original
     * MCFBlock: then, this MCFBlock may be copied from a MCFBlock that has
     * closed or deleted arcs and this method would not work. */

    /*void add_Modification(sp_Mod &mod) override
    {
      if (std::dynamic_pointer_cast<const NBModification>(mod))
      {
        // this is the "nuclear option": the MCFBlock has been re-loaded, so
        // the MCFClass solver also has to (immediately)
        //TODO: change MCFC function to Algo function.
        auto MCFB = static_cast<MCFBlock *>(f_Block);
        MCFC::LoadNet(MCFB->get_MaxNNodes(), MCFB->get_MaxNArcs(),
                      MCFB->get_NNodes(), MCFB->get_NArcs(),
                      MCFB->get_U().empty() ? nullptr : MCFB->get_U().data(),
                      MCFB->get_C().empty() ? nullptr : MCFB->get_C().data(),
                      MCFB->get_B().empty() ? nullptr : MCFB->get_B().data(),
                      MCFB->get_SN().data(), MCFB->get_EN().data());
        // TODO: PreProcess() changes the internal data of the MCFSolver using
        //       information about how the data of the MCF is *now*. If the
        //       data changes, some of the deductions (say, reducing the capacity
        //       of some arcs) may no longer be correct and they should be undone,
        //       but there isn't any proper way to handle this. Thus, PreProcess()
        //       has to be disabled for now; maybe later on someone will take
        //       care to make this work (or maybe not).
        // MCFC::PreProcess();
        // besides, any outstanding modification makes no sense any longer
        mod_clear();
      }
      else
        push_back(mod);
    }*/

    /*--------------------------------------------------------------------------*/
    /*-------------------------------- FRIENDS ---------------------------------*/
    /*--------------------------------------------------------------------------*/

    friend class MCFLemonState; // make MCFSolverState friend

    /** @} ---------------------------------------------------------------------*/
    /*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
    /*--------------------------------------------------------------------------*/

  protected:
    /*--------------------------------------------------------------------------*/
    /*-------------------------- PROTECTED METHODS -----------------------------*/
    /*--------------------------------------------------------------------------*/

    void process_outstanding_Modification(void);

    void guts_of_poM(c_p_Mod mod);

    /*--------------------------------------------------------------------------*/
    /*---------------------------- PROTECTED FIELDS  ---------------------------*/
    /*--------------------------------------------------------------------------*/

    const static std::vector<int> Solver_2_MCFClass_int;
    // the (static const) map between Solver int parameters and MCFClass ones

    const static std::vector<int> Solver_2_MCFClass_dbl;
    // the (static const) map between Solver int parameters and MCFClass ones

    std::string f_dmx_file; 
    // string for DMX file output


    /*--------------------------------------------------------------------------*/
    /*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
    /*--------------------------------------------------------------------------*/

  private:
    /*--------------------------------------------------------------------------*/
    /*-------------------------- PRIVATE METHODS -------------------------------*/
    /*--------------------------------------------------------------------------*/

    SMSpp_insert_in_factory_h;

     /*--------------------------------------------------------------------------*/
    /*-------------------------- PRIVATE FIELDS -------------------------------*/
    /*--------------------------------------------------------------------------*/
    int status = SMSpp_di_unipi_it::LEMON_sol_type::UNSOLVED;  //Variable used in compute function for getting status
    typename SMSppCostScaling< GR , V , C >::ProblemType status_2_pType;
    //Status of compute() method
    int f_method_type; //Variable used in set_par for switching value with f_method
    double ticks;  //Elapsed time in ticks for compute() method

    /*--------------------------------------------------------------------------*/

  }; // end( class MCFLemonSolver<SMSppCostScaling, GR, V, C>  Specialization)



/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register the various LEMONSolver< Alg , GR , V , C > to the Solver factory

  SMSpp_insert_in_factory_cpp_0_t(
 MCFLemonSolver< NetworkSimplex , SmartDigraph , double , double >);



/*--------------------------------------------------------------------------*/
/*----------------------- METHODS of MCFLemonSolver ------------------------*/
/*--------------------------------------------------------------------------*/

<<<<<<< HEAD
template <template <typename, typename, typename> class Algo, LEMONGraph GR,
          typename V, typename C>
void MCFLemonSolver<Algo, GR, V, C>::guts_of_set_Block(MCFBlock *MCFB) {
=======


/*--------------------------------------------------------------------------*/
/* Managing parameters for MCFSimplex---------------------------------------*/
/*
 * MCFSimplex has the following extra parameters:
 *
 * - kAlgPrimal     parameter to set algorithm (Primal/Dual):
 * - kAlgPricing    parameter to set algorithm of pricing
 * - kNumCandList   parameter to set the number of candidate list for
 *                  Candidate List Pivot method
 * - kHotListSize   parameter to set the size of Hot List
 *
 * (plus, actually, a few about Quadratic MCF that are not relevent here).
 * These are all "int" parameters, hence the "double" versions only issue the
 * method of the base CDASolver class, and therefore need not be defined. */
>>>>>>> 3587a0f95544d63e767d768c535e345b48386907

    // delete previous graph and f_algo (if present)
    // delete dgp;
    delete f_algo;
    f_algo = nullptr;
    // create and clear new Graph (MCFListDigraph or SmartDigraph)

    dgp = new GR;

    dgp->clear();

    // actually reserve get_MaxNNodes() node space in dgp
    auto maxNNodes = MCFB->get_MaxNNodes();
    dgp->reserveNode(maxNNodes);
    MCFBlock::Index n = MCFB->get_NNodes();
    // add all the nodes, one by one (looks stupid you con't do all in one blow,
    // but who are we to say ...)
    for (MCFBlock::Index i = 0; i < n; ++i)
        dgp->addNode();

    // reserve get_MaxNArcs() arcs space
    auto maxNArcs = MCFB->get_MaxNArcs();
    dgp->reserveArc(maxNArcs);
    MCFBlock::Index m = MCFB->get_NArcs();
    n_arcs = m;
    n_arcs_added = 0;
    // get starting node subset and ending node subset
    MCFBlock::c_Subset &sn = MCFB->get_SN();
    MCFBlock::c_Subset &en = MCFB->get_EN();

    // add all arcs in dgp
    // sn[ i ] - 1 and en[ i ] - 1 are used because sn , en nodes start from 1
    for (MCFBlock::Index i = 0; i < m; ++i)
        if constexpr (std::is_same<GR, SmartDigraph>::value) {
            dgp->addArc(dgp->nodeFromId(sn[i] - 1), dgp->nodeFromId(en[i] - 1));
        } else {
            dgp->addArc(sn[i] - 1, en[i] - 1);
        }

    // new instance of Lemon Algorithm with no arc mixing
    if constexpr (std::is_same<Algo<GR, V, C>,
                               NetworkSimplex<GR, V, C>>::value) {
        f_algo = new Algo<GR, V, C>(*dgp, false);
    } else {
        f_algo = new Algo<GR, V, C>(*dgp);
    }

    // Algo<GR, V, C> * algo = f_algo->reset();
    // defining names for types for readability
    using MCFArcMapV = typename GR::template ArcMap<V>;
    using MCFNodeMapV = typename GR::template NodeMap<V>;

    // now we are going to fill up ArcMap and NodeMap
    // this is the case of upperMap
    if (!MCFB->get_U().empty()) {
        um = new MCFArcMapV(*dgp); // create the upper map
        auto &u = MCFB->get_U();
        for (MCFBlock::Index i = 0; i < m; ++i)
            um->set(dgp->arcFromId(i), u[i]);

        f_algo->upperMap(*um); // pass it to the Algo
    }

    // this is the case of CostMap
    if (!MCFB->get_C().empty()) {
        cm = new MCFArcMapV(*dgp); // create the cost map
        auto &c = MCFB->get_C();
        for (MCFBlock::Index i = 0; i < m; ++i)
            cm->set(dgp->arcFromId(i), c[i]);

        f_algo->costMap(*cm); // pass it to the Algo
    }

    // this is the case of supplyMap
    if (!MCFB->get_B().empty()) {
        bm = new MCFNodeMapV(*dgp); // create the supply map
        auto &b = MCFB->get_B();
        for (MCFBlock::Index i = 0; i < n; i++)
            bm->set(dgp->nodeFromId(i),
                    b[i] == 0 ? 0 : -b[i]); // note tat supply = - deficit

        f_algo->supplyMap(*bm); // pass it to the algo
    }
}

template <template <typename, typename, typename> class Algo, LEMONGraph GR,
          typename V, typename C>
void MCFLemonSolver<Algo, GR, V, C>::set_Block(Block *block) {
    if (block == f_Block) // actually doing nothing
        return;           // cowardly and silently return

    Solver::set_Block(block); // attach to the new Block

    if (!block) // this is just: go sit down in a corner and wait
        return; // all done

    auto MCFB = dynamic_cast<MCFBlock *>(block);
    if (!MCFB)
        throw(std::invalid_argument(
            "MCFSolver:set_Block: block must be a MCFBlock"));

    bool owned = MCFB->is_owned_by(f_id);
    if ((!owned) && (!MCFB->read_lock()))
        throw(std::logic_error("cannot acquire read_lock on MCFBlock"));

    guts_of_set_Block(MCFB); // fill up the graph and f_algo

    if (!owned)
        MCFB->read_unlock();

} // end( MCFLemonSolver< Algo , GR , V , C >set_Block )

/*--------------------------------------------------------------------------*/

template <template <typename, typename, typename> class Algo, LEMONGraph GR,
          typename V, typename C>
int MCFLemonSolver<Algo, GR, V, C>::compute(bool changedvars) {

    using MCFArcMapV = typename GR::template ArcMap<V>;
    using MCFNodeMapV = typename GR::template NodeMap<V>;

    const static std::array<int, 3> LEMONstatus_2_sol_type = {
        Solver::kInfeasible, Solver::kOK, Solver::kUnbounded};

    lock(); // first of all, acquire self-lock

    if (!f_Block)              // there is no [MCFBlock] to solve
        return (kBlockLocked); // return error

    bool owned = f_Block->is_owned_by(f_id); // check if already locked
    if ((!owned) && (!f_Block->read_lock())) // if not try to read_lock
        return (kBlockLocked);               // return error on failure

    // while [read_]locked, process any outstanding Modification
    process_outstanding_Modification();

    if (!f_dmx_file.empty()) { // if so required
        // output the current instance (after the changes) to a DMX file
        std::ofstream ProbFile(f_dmx_file, ios_base::out | ios_base::trunc);
        if (!ProbFile.is_open())
            throw(std::logic_error("cannot open DMX file " + f_dmx_file));
        if constexpr (std::is_same<GR, SmartDigraph>::value) {
            writeDimacsMat(ProbFile, *dgp);
        }
        ProbFile.close();
    }

    if (!owned)                 // if the [MCF]Block was actually read_locked
        f_Block->read_unlock(); // read_unlock it

    // if the graph was modified adding\removing arcs, then reset() f_algo must
    // be used
    if (n_arcs_added > 0 || n_arcs_deleted > 0) {
        n_arcs += n_arcs_added;
        n_arcs -= n_arcs_deleted;
        n_arcs_added = 0;
        n_arcs_deleted = 0;

        f_algo->reset();

        cost_changed = true;
        cap_changed = true;
        supply_changed = true;
    }

    if (cost_changed) {
        f_algo->costMap(*cm);
        cost_changed = false;
    }

    if (cap_changed) {
        f_algo->upperMap(*um);
        cap_changed = false;
    }

    if (supply_changed) {
        f_algo->supplyMap(*bm);
        supply_changed = false;
    }

    // using Lemon Digraph Writer for debuggging
    // std::ofstream ProbFile("lmn.txt", ios_base::out | ios_base::trunc);
    // if (!ProbFile.is_open()) {
    //     throw(std::logic_error("cannot open file lmn.txt"));
    // }
    // DigraphWriter<GR>(*dgp, ProbFile).nodeMap("supply", *bm).arcMap("cost",
    // *cm).arcMap("upper", *um).run(); ProbFile.close();

    auto start = chrono::system_clock::now();
    guts_of_compute(); // here the actual magic is done by specialised classes
    auto end = chrono::system_clock::now();

    chrono::duration<double> elapsed = end - start;
    ticks = elapsed.count();

    unlock(); // release self-lock

    // now give out the result
    return (LEMONstatus_2_sol_type[this->get_status()]);

} // end( MCFLemonSolver< Algo , GR , V , C >::compute )

template <template <typename, typename, typename> class Algo, LEMONGraph GR,
          typename V, typename C>
void MCFLemonSolver<Algo, GR, V, C>::add_Modification(sp_Mod &mod) {

    if (std::dynamic_pointer_cast<const NBModification>(mod)) {
        // this is the "nuclear option": the MCFBlock has been re-loaded, so
        // the MCFClass solver also has to (immediately)
        auto MCFB = static_cast<MCFBlock *>(f_Block);
        guts_of_set_Block(MCFB);
        // besides, any outstanding modification makes no sense any longer
        mod_clear();
    } else
        push_back(mod);
}

template <template <typename, typename, typename> class Algo, LEMONGraph GR,
          typename V, typename C>
void MCFLemonSolver<Algo, GR, V, C>::process_outstanding_Modification(void) {

    // no-frills loop: do them in order, with no attempt at optimizing
    // note that NBModification have already been dealt with and therefore need
    // not be considered here

    for (;;) {
        auto mod = pop();
        if (!mod)
            break;

        guts_of_poM(mod.get());
    }
} // end( MCFLemonSolver::process_outstanding_Modification )

/*--------------------------------------------------------------------------*/

<<<<<<< HEAD
template <template <typename, typename, typename> class Algo, LEMONGraph GR,
          typename V, typename C>
void MCFLemonSolver<Algo, GR, V, C>::guts_of_poM(c_p_Mod mod) {
=======
template<>
Solver::idx_type MCFSolver< MCFSimplex >::get_num_int_par( void ) const
{
 //return( CDASolver::get_num_int_par() + 5 );
 }
>>>>>>> 3587a0f95544d63e767d768c535e345b48386907

    auto MCFB = static_cast<MCFBlock *>(f_Block);

    // process Modification - - - - - - - - - - - - - - - - - - - - - - - - - -
    // -
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //-
    /* This requires to patiently sift through the possible Modification types
     * to find what this Modification exactly is, and call the appropriate
     * method of LEMON. */

    // GroupModification- - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // -
    if (auto tmod = dynamic_cast<const GroupModification *>(mod)) {
        for (const auto &submod : tmod->sub_Modifications())
            guts_of_poM(submod.get());

<<<<<<< HEAD
        return;
    }

    // MCFBlockRngdMod- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // -
    /* Note: in the following we can assume that C, B and U are nonempty. This
     * is because they can be empty only if they are so when the object is
     * loaded. But if a Modification has been issued they are no longer empty (a
     * Modification changin nothing from the "empty" state is not issued). */
=======
template<>
int MCFSolver< MCFSimplex >::get_dflt_int_par( idx_type par ) const
{/*
 static const std::array< int , 5 > my_dflt_int_par = { MCFClass::kYes ,
		  MCFClass::kYes , MCFSimplex::kCandidateListPivot , 0 , 0 };

 return( par >= intLastParCDAS ? my_dflt_int_par[ par - intLastParCDAS ]
	                       : CDASolver::get_dflt_int_par( par ) );
 */}
>>>>>>> 3587a0f95544d63e767d768c535e345b48386907

    if (auto tmod = dynamic_cast<const MCFBlockRngdMod *>(mod)) {
        auto rng = tmod->rng();

        switch (tmod->type()) {
        case (MCFBlockMod::eChgCost): {

            if (rng.second == rng.first + 1) {

<<<<<<< HEAD
                if (dgp->valid(dgp->arcFromId(rng.first))) {
                    auto cost = MCFB->get_C(rng.first);
                    cm->set(dgp->arcFromId(rng.first),
                            std::isnan(cost) ? 0 : cost);
                    cost_changed = true;
                } else {
                    for (MCFBlock::Index i = rng.first; i < rng.second; i++) {
                        if (dgp->valid(dgp->arcFromId(i)) &&
                            !(std::isnan(MCFB->get_C(i)))) {
                            auto cost = MCFB->get_C(i);
                            cm->set(dgp->arcFromId(i), cost);
                            cost_changed = true;
                        } else if (std::isnan(MCFB->get_C(i))) {
                            cm->set(dgp->arcFromId(i), 0);
                            cost_changed = true;
                        }
                    }
                }

            } else {
                for (MCFBlock::Index i = rng.first; i < rng.second; i++) {
                    if (dgp->valid(dgp->arcFromId(i))) {
                        auto cost = MCFB->get_C(i);
                        cm->set(dgp->arcFromId(i), std::isnan(cost) ? 0 : cost);
                        cost_changed = true;
                    }
                }
            }
=======
template<>
Solver::idx_type MCFSolver< MCFSimplex >::int_par_str2idx(
					     const std::string & name ) const
{/*
 if( name == "kReopt" )
  return( intLastParCDAS );
 if( name == "kAlgPrimal" )
  return( intLastParCDAS + 1 );
 if( name == "kAlgPricing" )
  return( intLastParCDAS + 2 );
 if( name == "kNumCandList" )
  return( intLastParCDAS + 3 );
 if( name == "kHotListSize" )
  return( intLastParCDAS + 4 );

 return( CDASolver::dbl_par_str2idx( name ) );
 */}
>>>>>>> 3587a0f95544d63e767d768c535e345b48386907

            return;
        }

        /*--------------------------------------------------------------*/
        /*--------------------------------------------------------------*/
        case (MCFBlockMod::eChgCaps): {

            if (rng.second == rng.first + 1) {
                if (dgp->valid(dgp->arcFromId(rng.first))) {
                    um->set(dgp->arcFromId(rng.first), MCFB->get_U(rng.first));
                    cap_changed = true;
                }
            } else {
                for (MCFBlock::Index i = rng.first; i < rng.second; i++) {
                    if (dgp->valid(dgp->arcFromId(i))) {
                        um->set(dgp->arcFromId(i), MCFB->get_U(i));
                        cap_changed = true;
                    }
                }
            }

<<<<<<< HEAD
            return;
        }

        /*--------------------------------------------------------------*/
        /*--------------------------------------------------------------*/
        case (MCFBlockMod::eChgDfct): {
=======
template<>
const std::string & MCFSolver< MCFSimplex >::int_par_idx2str( idx_type idx )
 const {/*
 static const std::array< std::string , 5 > my_int_pars_str = {
  "kReopt" , "kAlgPrimal" , "kAlgPricing" , "kNumCandList" , "kHotListSize"
  };

 return( idx >= intLastParCDAS ? my_int_pars_str[ idx - intLastParCDAS ]
	                       : CDASolver::int_par_idx2str( idx ) );
 */}
>>>>>>> 3587a0f95544d63e767d768c535e345b48386907

            if (rng.second == rng.first + 1) {
                bm->set(dgp->nodeFromId(rng.first), -MCFB->get_B(rng.first));
                supply_changed = true;

            }

            else {
                for (MCFBlock::Index i = rng.first; i < rng.second; i++) {
                    bm->set(dgp->nodeFromId(i), -MCFB->get_B(i));
                    supply_changed = true;
                }
            }
            return;
        }

        /*--------------------------------------------------------------*/
        /*--------------------------------------------------------------*/
        case (MCFBlockMod::eOpenArc): {
            // only if the graph is ListDigraph
            if constexpr (std::is_same<GR, MCFListDigraph>::value) {

                for (; rng.first < rng.second; ++rng.first) {
                    if ((!MCFB->is_deleted(rng.first)) &&
                        dgp->valid(dgp->arcFromId(rng.first)) &&
                        dgp->isClosed(rng.first)) {
                        dgp->openArc(rng.first);
                        n_arcs_added++;
                    }
                }
            } else {
                throw(std::logic_error(
                    "SmartDigraph doesn't support operations on arc"));
            }
            return;
        }

        /*--------------------------------------------------------------*/
        /*--------------------------------------------------------------*/
        case (MCFBlockMod::eCloseArc): {
            // only if the graph is ListDigraph
            if constexpr (std::is_same<GR, MCFListDigraph>::value) {

                for (; rng.first < rng.second; ++rng.first) {
                    if ((!MCFB->is_deleted(rng.first)) &&
                        !dgp->isClosed(rng.first) &&
                        dgp->valid(dgp->arcFromId(rng.first))) {
                        dgp->closeArc(rng.first);
                        n_arcs_deleted++;
                    }
                }
            } else {
                throw(std::logic_error(
                    "SmartDigraph doesn't support operations on arc"));
            }
            return;
        }

        /*--------------------------------------------------------------*/
        /*--------------------------------------------------------------*/
        case (MCFBlockMod::eAddArc): {
            // only if the graph is ListDigraph
            if constexpr (std::is_same<GR, MCFListDigraph>::value) {

                auto ca = MCFB->get_C(rng.first);
                auto startNode = MCFB->get_SN(rng.first) - 1;
                auto endNode = MCFB->get_EN(rng.first) - 1;
                auto capacity = MCFB->get_U(rng.first);

                dgp->addArc(startNode, endNode);
                auto arc = dgp->arcFromId(rng.first);

                n_arcs_added++;
                cm->set(arc, std::isnan(ca) ? 0 : ca);
                um->set(arc, capacity);

                cost_changed = true;
                cap_changed = true;

                return;
            } else {
                throw(std::logic_error(
                    "SmartDigraph doesn't support operations on arc"));
            }

            return;
        }
        /*--------------------------------------------------------------*/
        /*--------------------------------------------------------------*/
        case (MCFBlockMod::eRmvArc): {
            // only if ListDigraph is used
            if constexpr (std::is_same<GR, MCFListDigraph>::value) {
                // if the arc was previous deleted, return
                if (!dgp->valid(dgp->arcFromId(rng.second - 1))) {
                    return;
                }
                auto arc = dgp->arcFromId(rng.second - 1);

                // if the arc isClosed, only 'mark' it as eliminated, otherwise
                // normally erase it
                if (dgp->isClosed(rng.second - 1)) {
                    dgp->eraseClosed(rng.second - 1);
                } else {
                    dgp->erase(arc);
                }

                n_arcs_deleted++;
            } else {
                throw(std::logic_error(
                    "SmartDigraph doesn't support operations on arc"));
            }

            return;
        }
        default:
            throw(std::invalid_argument("unknown MCFBlockRngdMod type"));
        }
    }

    // MCFBlockSbstMod- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // -
    if (auto tmod = dynamic_cast<const MCFBlockSbstMod *>(mod)) {
        switch (tmod->type()) {

        case (MCFBlockMod::eOpenArc): {
            // only if ListDigraph is used
            if constexpr (std::is_same<GR, MCFListDigraph>::value) {
                for (auto arc : tmod->nms())
                    // checking if the arc is closed
                    if ((!MCFB->is_deleted(arc)) && dgp->isClosed(arc) &&
                        dgp->valid(dgp->arcFromId(arc))) {
                        dgp->openArc(arc);
                        n_arcs_added++;
                    }
            } else {
                throw(std::logic_error(
                    "SmartDigraph doesn't support operations on arc"));
            }
            return;
        }

        /*--------------------------------------------------------------*/
        /*--------------------------------------------------------------*/
        case (MCFBlockMod::eCloseArc): {
            // only if ListDigraph is used
            if constexpr (std::is_same<GR, MCFListDigraph>::value) {

                for (auto arc : tmod->nms())
                    // checking if the arc is not closed
                    if ((!MCFB->is_deleted(arc)) &&
                        dgp->valid(dgp->arcFromId(arc)) &&
                        !dgp->isClosed(arc)) {
                        dgp->closeArc(arc);
                        n_arcs_deleted++;
                    }

                return;
            } else {
                throw(std::logic_error(
                    "SmartDigraph doesn't support operations on arc"));
            }
        };

        /*--------------------------------------------------------------*/
        /*--------------------------------------------------------------*/
        case (MCFBlockMod::eChgCost): {

            for (auto i : tmod->nms())
                if (dgp->valid(dgp->arcFromId(i))) {
                    auto arc = dgp->arcFromId(i);
                    auto cost = MCFB->get_C(i);
                    cm->set(arc, std::isnan(cost) ? 0 : cost);
                    cost_changed = true;
                }

            return;
        }

        /*--------------------------------------------------------------*/
        /*--------------------------------------------------------------*/
        case (MCFBlockMod::eChgCaps): {

            auto &CC = MCFB->get_C();
            auto &U = MCFB->get_U();
            for (auto i : tmod->nms())
                if (!std::isnan(CC[i]) && dgp->valid(dgp->arcFromId(i))) {

                    um->set(dgp->arcFromId(i), U[i]);
                    cap_changed = true;
                }

            return;
        }

        /*--------------------------------------------------------------*/
        /*--------------------------------------------------------------*/
        case (MCFBlockMod::eChgDfct): {
            MCFBlock::Vec_FNumber NDfct(tmod->nms().size());
            MCFBlock::Subset nmsI(tmod->nms().size() + 1);
            *copy(tmod->nms().begin(), tmod->nms().end(), nmsI.begin()) =
                Inf<MCFBlock::Index>();
            auto B = MCFB->get_B();
            for (MCFBlock::Index i = 0; i < NDfct.size(); i++) {
                bm->set(dgp->nodeFromId(i), -B[nmsI[i]]);
                supply_changed = true;
            }
            return;
        }

        default:
            throw(std::invalid_argument("unknown MCFBlockSbstMod type"));
        }
    }

    // any remaining Modification is plainly ignored, since it must be an
    // "abstract" Modification, which this Solver does not need to look at

} // end( guts_of_poM )

/*--------------------------------------------------------------------------*/
<<<<<<< HEAD
/*------------------- End File MCFLemonSolver.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
=======
/* Managing parameters for RelaxIV -----------------------------------------*/
/*
 * RelaxIV has the following extra parameters:
 *
 * - kAuction     the auction/shortest paths initialization is used
 *
 * These are all "int" parameters, hence the "double" versions only issue the
 * method of the base CDASolver class, and therefore need not be defined. */

#ifdef HAVE_RELAX

/*--------------------------------------------------------------------------*/

template<>
const std::vector< int > MCFSolver< RelaxIV >::Solver_2_MCFClass_int = {
 MCFClass::kMaxIter ,        // intMaxIter
 -1 ,                        // intMaxThread
 -1 ,                        // intEverykIt
 -1 ,                        // intMaxSol
 -1 ,                        // intLogVerb
 -1 ,                        // intMaxDSol
 MCFClass::kReopt ,          // intLastParCDAS
 RelaxIV::kAuction
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<>
const std::vector< int > MCFSolver< RelaxIV >::Solver_2_MCFClass_dbl = {
 MCFClass::kMaxTime ,       // dblMaxTime
 -1 ,                       // dblEveryTTm
 -1 ,                       // dblRelAcc
 MCFClass::kEpsFlw ,        // dblAbsAcc
 -1 ,                       // dblUpCutOff
 -1 ,                       // dblLwCutOff
 -1 ,                       // dblRAccSol
 -1 ,                       // dblAAccSol
 -1 ,                       // dblFAccSol
 -1 ,                       // dblRAccDSol
 MCFClass::kEpsCst ,        // dblAAccDSol
 -1                         // dblFAccDSol
 };

/*--------------------------------------------------------------------------*/

template<>
Solver::idx_type MCFSolver< RelaxIV >::get_num_int_par( void ) const
{
 //return( CDASolver::get_num_int_par() + 2 );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
Solver::idx_type MCFSolver< RelaxIV >::get_num_dbl_par( void ) const { }

----------------------------------------------------------------------------*/

template<>
int MCFSolver< RelaxIV >::get_dflt_int_par( const idx_type par ) const
{/*
 static const std::array< int , 2 > my_dflt_int_par = { MCFClass::kYes ,
						        MCFClass::kYes };

 return( par >= intLastParCDAS ? my_dflt_int_par[ par - intLastParCDAS ]
	                       : CDASolver::get_dflt_int_par( par ) );
 */}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
double MCFSolver< RelaxIV >::get_dflt_dbl_par( idx_type par ) const { }

----------------------------------------------------------------------------*/

template<>
Solver::idx_type MCFSolver< RelaxIV >::int_par_str2idx(
					     const std::string & name ) const
{/*
 if( name == "kReopt" )
  return( intLastParCDAS );
 if( name == "kAuction" )
  return( intLastParCDAS + 1 );

 return( CDASolver::dbl_par_str2idx( name ) );
 */}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
Solver::idx_type MCFSolver< RelaxIV >::dbl_par_str2idx(
					 const std::string & name ) const { }

----------------------------------------------------------------------------*/

template<>
const std::string & MCFSolver< RelaxIV >::int_par_idx2str( idx_type idx )
 const {
 /*static const std::array< std::string , 2 > my_int_pars_str = { "kReopt" ,
								"kAuction" };

 return( idx >= intLastParCDAS ? my_int_pars_str[ idx - intLastParCDAS ]
	                       : CDASolver::int_par_idx2str( idx ) );
 */}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
const std::string & MCFSolver< RelaxIV >::dbl_par_idx2str( idx_type idx )
 const { }
 */

#endif

/*--------------------------------------------------------------------------*/
/* Managing parameters for CPLEX -------------------------------------------*/

#ifdef HAVE_CPLEX

/*--------------------------------------------------------------------------*/

template<>
const std::vector< int > MCFSolver< MCFCplex >::Solver_2_MCFClass_int = {
 MCFClass::kMaxIter ,        // intMaxIter
 -1 ,                        // intMaxThread
 -1 ,                        // intEverykIt
 -1 ,                        // intMaxSol
 -1 ,                        // intLogVerb
 -1 ,                        // intMaxDSol
 MCFClass::kReopt ,          // intLastParCDAS
 MCFCplex::kQPMethod
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<>
const std::vector< int > MCFSolver< MCFCplex >::Solver_2_MCFClass_dbl =
{
 MCFClass::kMaxTime ,       // dblMaxTime
 -1 ,                       // dblEveryTTm
 -1 ,                       // dblRelAcc
 MCFClass::kEpsFlw ,        // dblAbsAcc
 -1 ,                       // dblUpCutOff
 -1 ,                       // dblLwCutOff
 -1 ,                       // dblRAccSol
 -1 ,                       // dblAAccSol
 -1 ,                       // dblFAccSol
 -1 ,                       // dblRAccDSol
 MCFClass::kEpsCst ,        // dblAAccDSol
 -1                         // dblFAccDSol
 };

/*--------------------------------------------------------------------------*/

template<>
Solver::idx_type MCFSolver< MCFCplex >::get_num_int_par( void ) const {
 //return( CDASolver::get_num_int_par() + 2 );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
Solver::idx_type MCFSolver< MCFCplex >::get_num_dbl_par( void ) const { }

----------------------------------------------------------------------------*/

template<>
int MCFSolver< MCFCplex >::get_dflt_int_par( idx_type par ) const {
 /*static const std::array< int , 2 > my_dflt_int_par = { MCFClass::kYes ,
                                                        MCFClass::kYes };
 return( par >= intLastParCDAS ?
         my_dflt_int_par[ par - intLastParCDAS ] :
         CDASolver::get_dflt_int_par( par ) );
 */}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 template<>
 double MCFSolver< MCFCplex >::get_dflt_dbl_par( idx_type par ) const { }

----------------------------------------------------------------------------*/

template<>
Solver::idx_type MCFSolver< MCFCplex >::int_par_str2idx(
					   const std::string & name ) const {
 /*if( name == "kReopt" )
  return( intLastParCDAS );
 if( name == "kQPMethod" )
  return( kQPMethod + 1 );

 return( CDASolver::dbl_par_str2idx( name ) );
 */}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
Solver::idx_type MCFSolver< MCFCplex >::dbl_par_str2idx(
                                         const std::string & name ) const { }

----------------------------------------------------------------------------*/

template<>
const std::string & MCFSolver< MCFCplex >::int_par_idx2str( idx_type idx )
 const {
 /*static const std::array< std::string , 2 > my_int_pars_str = { "kReopt" ,
                                                                "kQPMethod" };
 return( idx >= intLastParCDAS ?
         my_int_pars_str[ idx - intLastParCDAS ] :
         CDASolver::int_par_idx2str( idx ) );
 */}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
const std::string & MCFSolver< MCFCplex >::dbl_par_idx2str( idx_type idx )
 const { }
*/

#endif

/*--------------------------------------------------------------------------*/
/* Managing parameters for SPTree ------------------------------------------*/

#ifdef HAVE_SPTRE

/*--------------------------------------------------------------------------*/

template<>
const std::vector< int > MCFSolver< SPTree >::Solver_2_MCFClass_int = {
 MCFClass::kMaxIter,        // intMaxIter
 -1,                        // intMaxSol
 -1,                        // intLogVerb
 -1,                        // intMaxDSol
 MCFClass::kReopt,          // intLastParCDAS
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<>
const std::vector< int > MCFSolver< SPTree >::Solver_2_MCFClass_dbl = {
 MCFClass::kMaxTime,        // dblMaxTime
 -1,                        // dblRelAcc
 MCFClass::kEpsFlw,         // dblAbsAcc
 -1,                        // dblUpCutOff
 -1,                        // dblLwCutOff
 -1,                        // dblRAccSol
 -1,                        // dblAAccSol
 -1,                        // dblFAccSol
 -1,                        // dblRAccDSol
 MCFClass::kEpsCst,         // dblAAccDSol
 -1                         // dblFAccDSol
 };

/*--------------------------------------------------------------------------*/

template<>
Solver::idx_type MCFSolver< SPTree >::get_num_int_par( void ) const {
 /*return( CDASolver::get_num_int_par() + 1 );
 */}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
Solver::idx_type MCFSolver< SPTree >::get_num_dbl_par( void ) const { }

----------------------------------------------------------------------------*/

template<>
int MCFSolver< SPTree >::get_dflt_int_par( idx_type par ) const {
 /*static const std::array< int , 1 > my_dflt_int_par = { MCFClass::kYes };

 return( par >= intLastParCDAS ?
         my_dflt_int_par[ par - intLastParCDAS ] :
         CDASolver::get_dflt_int_par( par ) );
 */}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
double MCFSolver< SPTree >::get_dflt_dbl_par( idx_type par ) const { }

----------------------------------------------------------------------------*/

template<>
Solver::idx_type MCFSolver< SPTree >::int_par_str2idx(
					  const std::string & name ) const {
 /*if( name == "kReopt" )
  return( intLastParCDAS );

 return( CDASolver::dbl_par_str2idx( name ) );
 */}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
Solver::idx_type MCFSolver< SPTree >::dbl_par_str2idx(
                                         const std::string & name ) const { }

----------------------------------------------------------------------------*/

template<>
const std::string & MCFSolver< SPTree >::int_par_idx2str( idx_type idx )
 const {/*
 static const std::array< std::string , 1 > my_int_pars_str = { "kReopt" };

 return( idx >= intLastParCDAS ?
         my_int_pars_str[ idx - intLastParCDAS ] :
         CDASolver::int_par_idx2str( idx ) );
 */}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 template<>
 const std::string & MCFSolver< SPTree >::dbl_par_idx2str( idx_type idx )
 const { }
*/

#endif

/*--------------------------------------------------------------------------*/
/*----------------------- End File MCFSolver.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
 
>>>>>>> 3587a0f95544d63e767d768c535e345b48386907

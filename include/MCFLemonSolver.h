/*--------------------------------------------------------------------------*/
/*------------------------ File MCFLemonSolver.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the MCFLemonSolver class, implementing the Solver
 * interface, in particular in its CDASolver version, for Min-Cost Flow
 * problems as set by MCFBlock.
 *
 * This is based on interfacing algorithms implemented in the LEMON
 * (Library for Efficient Modeling and Optimization in Networks) project,
 * as currently found at
 *
 *     https://lemon.cs.elte.hu/trac/lemon
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
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __MCFLemonSolver
<<<<<<< HEAD
#define __MCFLemonSolver
/* self-identification: #endif at the end of the file */
=======
 #define __MCFLemonSolver
                      /* self-identification: #endif at the end of the file */
>>>>>>> 3587a0f95544d63e767d768c535e345b48386907

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "CDASolver.h"

#include "MCFBlock.h"

#include "lemon/capacity_scaling.h"

#include "lemon/cost_scaling.h"

#include "lemon/cycle_canceling.h"

#include "lemon/network_simplex.h"

//!!
#include <concepts>

#include <lemon/list_graph.h>

#include <lemon/smart_graph.h>

#include <lemon/dimacs.h>

#include <lemon/concepts/digraph.h>

<<<<<<< HEAD
#include <type_traits>

#include <lemon/lgf_writer.h>

#include <utility>
=======
>>>>>>> 3587a0f95544d63e767d768c535e345b48386907

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/
/// namespace for the Structured Modeling System++ (SMS++)
<<<<<<< HEAD
namespace SMSpp_di_unipi_it {
using namespace lemon;
using namespace lemon::concepts;
using namespace std;

/** @} ---------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup MCFListDigraph_CLASSES Classes in MCFLemonSolver.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS MCFListDigraph ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/**
 * MCFListDigraph implements ListDigraph for customizing functions
 * that provides access to the structure of the graph used by f_algo. This class
 * comes from the needs of implementing openArc, closeArc method, useful in
 * Modification for the ListDigraph directed graph class provided by LEMON.
 * Keeping the arc "closed" while f_algo runs, it is possible only
 * changing to negative value the startID of the arc considered, while for
 * "open" an arc it is possible restoring the original value of the startID. An
 * implementation of eraseClosed is provided, since close() method modifies the
 * structure of the arc object without 'mark' it deleted, and since an arc can
 * be removed while it is closed, we need to 'mark' it with the standard LEMON
 * operations to recognise deleted arcs( that is, set prev_in field of n-th
 * arc to -2).
 *
 *
 * MCFArc is provided since constructor of the Arc class must be used in addArc
 * method, because it sets the value of the arcID in the arc object, required
 * for iterate over the graph structure of LEMON classes.
 */

class MCFListDigraph : public lemon::ListDigraph {
  public:
    MCFListDigraph() : lemon::ListDigraph() {}
    ~MCFListDigraph() = default;

    class MCFArc : ListDigraph::Arc {
        friend class MCFListDigraph;

      protected:
        MCFArc(int pid) { id = pid; }
    };

    using ListDigraphBase::Arc;
    using ListDigraphBase::first_free_arc;
    using ListDigraphBase::Node;

    bool isClosed(int n) { return _arcs[n].source < 0; }

    MCFArc addArc(int u, int v) {
        int n;

        if (first_free_arc == -1) {
            n = _arcs.size();
            _arcs.push_back(ArcT());
        } else {
            // look for the min-index arc in the deleted arc list
            int minidx = first_free_arc;
            int idx = first_free_arc;
            int prevmin = -1;
            for (auto nextidx = _arcs[idx].next_in; nextidx != -1; idx = nextidx, nextidx = _arcs[idx].next_in)
                if (nextidx < minidx) {
                    minidx = nextidx;
                    prevmin = idx;
                }

            if (prevmin == -1)
                first_free_arc = _arcs[first_free_arc].next_in;
            else
                _arcs[prevmin].next_in = _arcs[minidx].next_in;

            _arcs[minidx].source = u;
            _arcs[minidx].target = v;

            _arcs[minidx].next_out = _nodes[u].first_out;
            if (_nodes[u].first_out != -1) {
                _arcs[_nodes[u].first_out].prev_out = minidx;
            }

            _arcs[minidx].next_in = _nodes[v].first_in;
            if (_nodes[v].first_in != -1) {
                _arcs[_nodes[v].first_in].prev_in = minidx;
            }

            _arcs[minidx].prev_in = _arcs[minidx].prev_out = -1;

            _nodes[u].first_out = _nodes[v].first_in = minidx;

            return MCFArc(minidx);
        }
        // otherwise add to the last position
        _arcs[n].source = u;
        _arcs[n].target = v;

        _arcs[n].next_out = _nodes[u].first_out;
        if (_nodes[u].first_out != -1) {
            _arcs[_nodes[u].first_out].prev_out = n;
        }

        _arcs[n].next_in = _nodes[v].first_in;
        if (_nodes[v].first_in != -1) {
            _arcs[_nodes[v].first_in].prev_in = n;
        }

        _arcs[n].prev_in = _arcs[n].prev_out = -1;

        _nodes[u].first_out = _nodes[v].first_in = n;

        return MCFArc(n);
    }
    void closeArc(int n) {

        if (_arcs[n].next_in != -1) {
            _arcs[_arcs[n].next_in].prev_in = _arcs[n].prev_in;
        }

        if (_arcs[n].prev_in != -1) {
            _arcs[_arcs[n].prev_in].next_in = _arcs[n].next_in;
        } else {
            _nodes[_arcs[n].target].first_in = _arcs[n].next_in;
        }

        if (_arcs[n].next_out != -1) {
            _arcs[_arcs[n].next_out].prev_out = _arcs[n].prev_out;
        }

        if (_arcs[n].prev_out != -1) {
            _arcs[_arcs[n].prev_out].next_out = _arcs[n].next_out;
        } else {
            _nodes[_arcs[n].source].first_out = _arcs[n].next_out;
        }

        _arcs[n].source = -(_arcs[n].source + 1);
    }

    void openArc(int n) {

        _arcs[n].source = -(_arcs[n].source + 1);
        auto uid = _arcs[n].source;
        auto vid = _arcs[n].target;

        _arcs[n].next_out = _nodes[uid].first_out;
        if (_nodes[uid].first_out != -1) {
            _arcs[_nodes[uid].first_out].prev_out = n;
        }

        _arcs[n].next_in = _nodes[vid].first_in;
        if (_nodes[vid].first_in != -1) {
            _arcs[_nodes[vid].first_in].prev_in = n;
        }

        _arcs[n].prev_in = _arcs[n].prev_out = -1;

        _nodes[uid].first_out = _nodes[vid].first_in = n;
    }

    void eraseClosed(int n) {

        _arcs[n].next_in = first_free_arc;
        first_free_arc = n;
        _arcs[n].prev_in = -2;
    }
};
=======
namespace SMSpp_di_unipi_it
{
 using namespace lemon;
 using namespace lemon::concepts;
 using namespace std;
>>>>>>> 3587a0f95544d63e767d768c535e345b48386907

/*--------------------------------------------------------------------------*/
/*-------------------------- TEMPLATE TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup MCFLemonSolver_TYPES Template Types in MCFLemonSolver.h
 *
 * Algorithms in the LEMON projects are template over at least three types:
 *
 * - GR, which is the type of graph
<<<<<<< HEAD
 *
 *   The possibilities are:
 *
=======
 *  
 *   The possibilities are:
 * 
>>>>>>> 3587a0f95544d63e767d768c535e345b48386907
 *   - SmartDigraph is a simple and fast directed graph implementation
 *     It is quite memory efficient but at the price that it does not support
 *     node and arc deletion.
 *     It will be hard find the best way for modification.
<<<<<<< HEAD
 *
=======
 * 
>>>>>>> 3587a0f95544d63e767d768c535e345b48386907
 *   - ListDigraph, a versatile and fast directed graph implementation based
 *     on linked lists that are stored in std::vector structures.
 *     This class provides only linear time counting for nodes and arcs.
 *     It support node and arc deletion, useful for modification.
 *
 * - V, which is the type of flows / deficits; typically, double can be used
 *   for maximum compatibility, but int (or even smaller) would yeld better
 *   performances;
 *
 * - C, which is the type of ar costs; typically, double can be used for
 *   maximum compatibility, but int (or even smaller) would yeld better
 *   performances;
 *
 * Furthermore, scaling-type algorithms may behave in different ways
 * according to which combination of V and C is used, and there are different
 * "traits" for this which are another scaling parameter. However, we prefer
 * that the MCFLemonSolver class is always template over the three first
 * parameter only, which is why we define SMSppCapacityScaling and
 * SMSppCostScaling as template over < GR , V , C > and using the default
 * trait.
 *
 * Thus, the concept LEMONGraph is defined that only allows all possible
 * types of LEMON graphs, i.e.,
<<<<<<< HEAD
 *
=======
 * 
>>>>>>> 3587a0f95544d63e767d768c535e345b48386907
 *   - SmartDigraph is a simple and fast directed graph implementation
 *     It is quite memory efficient but at the price that it does not support
 *     node and arc deletion.
 *     It will be hard find the best way for modification.
<<<<<<< HEAD
 *
=======
 * 
>>>>>>> 3587a0f95544d63e767d768c535e345b48386907
 *   - ListDigraph, a versatile and fast directed graph implementation based
 *     on linked lists that are stored in std::vector structures.
 *     This class provides only linear time counting for nodes and arcs.
 *     It support node and arc deletion, useful for modification.
 *
 * Similarly, the concept LEMONAlgorithm is defined that only allows all
 * possible types of LEMON algorithms, i.e., [SMSpp]CapacityScaling,
 * [SMSpp]CostScaling, CycleCanceling, and NetworkSimplex.
 *
 *  @{ */

/// concept for "one of the LEMON graphs"
template <typename Type>
concept LEMONGraph = std::is_same<Type, SmartDigraph>::value || std::is_same<Type, MCFListDigraph>::value;
/// CapacityScaling algorithm using the default trait
template <LEMONGraph GR, typename V, typename C>
class SMSppCapacityScaling : public CapacityScaling<GR, V, C, CapacityScalingDefaultTraits<GR, V, C>> {
  public:
    SMSppCapacityScaling(const GR &dgp) : CapacityScaling<GR, V, C, CapacityScalingDefaultTraits<GR, V, C>>(dgp) {}

    ~SMSppCapacityScaling() = default;
};

/// CostScaling algorithm using the default trait
template <LEMONGraph GR, typename V, typename C>
class SMSppCostScaling : public CostScaling<GR, V, C, CostScalingDefaultTraits<GR, V, C>> {
  public:
    SMSppCostScaling(const GR &dgp) : CostScaling<GR, V, C, CostScalingDefaultTraits<GR, V, C>>(dgp) {}

<<<<<<< HEAD
    ~SMSppCostScaling() = default;
};

/** @} ---------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup MCFLemonSolver_CLASSES Classes in MCFLemonSolver.h
 *  @{ */
=======
 /// CostScaling algorithm using the default trait
 template < LEMONGraph GR , typename V , typename C>
 class SMSppCostScaling : public
  CostScaling< GR , V , C , CostScalingDefaultTraits< GR , V , C > >
  {};


>>>>>>> 3587a0f95544d63e767d768c535e345b48386907

/** @} ---------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
<<<<<<< HEAD
/*------------------------ CLASS MCFLemonSolver ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// CDASolver for MCFBlock based on the LEMON project
/** The MCFLemonSolver implements Solver interface for MCFBlock that represent
 * (Linear) Min-Cost Flow (MCF) problems, using algorithms of LEMON library.
 * Because MCF is a Linear Program it has a(n exact) dual, and therefore
 * MCFLemonSolver implements the CDASolver interface for also giving out dual
 * information.
 *
 * The MCFLemonSolver is template over four different types:
 *
 * - GR, which is the type of graph:
 *
 *   The possibilities are:
 *   - SmartDigraph is a simple and fast directed graph implementation
 *     It is quite memory efficient but at the price that it does not support
 *     node and arc deletion.
 *     It only supports Modification for costs, supplies and capacities.
 *
 *   - ListDigraph, a versatile and fast directed graph implementation based
 *     on linked lists that are stored in std::vector structures.
 *     This class provides only linear time counting for nodes and arcs.
 *     It support also node and arc deletion, useful for modification.
 *
 * - V, which is the type of flows / deficits; typically, double can be used
 *   for maximum compatibility, but int (or even smaller) would yeld better
 *   performances;
 *
 * - C, which is the type of ar costs; typically, double can be used for
 *   maximum compatibility, but int (or even smaller) would yeld better
 *   performances;
 *
 * - Algo, which is the specific algorithm (itself, template over GR, V, and
 *   C) implemented in the LEMON package. The possibilities are:
 *
 *   = NetworkSimplex implements the primal Network Simplex algorithm for
 *     finding a minimum cost flow. This algorithm is a highly efficient
 *     specialized version of the linear programming simplex method directly
 *     for the minimum cost flow problem.
 *
 *   = CycleCanceling implements three different cycle-canceling algorithms
 *     for finding a minimum cost flow. The most efficent one is the
 *     Cancel-and-tighten algorithm, thus it is the default method. It runs
 *     in strongly polynomial time, but in practice, it is typically orders of
 *     magnitude slower than the scaling algorithms and NetworkSimplex.
 *
 *   = CostScaling implements a cost scaling algorithm that performs
 *     push/augment and relabel operations for finding a minimum cost flow.
 *     It is a highly efficient primal-dual solution method, which can be
 *     viewed as the generalization of the preflow push-relabel algorithm for
 *     the maximum flow problem. It is a polynomial algorithm.
 *
 *   = CapacityScaling implements the capacity scaling version of the
 *     successive shortest path algorithm for finding a minimum cost flow. It
 *     is an efficient dual solution method, which runs in polynomial time.
 *     In special cases it can be more efficient than CostScaling and
 *     NetworkSimplex algorithms.
 *
 *   In general, NetworkSimplex and CostScaling are the fastest
 *   implementations available in LEMON for solving this problem.
 *
 *   Note that scaling-type algorithms may behave in different ways according
 *   to which combination of V and C is used, and there are different
 *   "traits" for this which are another scaling parameter. However, in order
 *   to make MCFLemonSolver class template over always the same number of
 *   template parameters we fix the use of the default trait, which is why
 *   SMSppCapacityScaling and SMSppCostScaling are defined (as template over
 *   < GR , V , C >) that are meant to be used instead of the original
 *   CapacityScaling and CostScaling. */

template <template <typename, typename, typename> class Algo, LEMONGraph GR, typename V, typename C>
class MCFLemonSolver : public CDASolver {
    /*--------------------------------------------------------------------------*/
    /*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
    /*--------------------------------------------------------------------------*/

  public:
    /*--------------------------------------------------------------------------*/
    /*---------------------------- PUBLIC TYPES --------------------------------*/
    /*--------------------------------------------------------------------------*/

    using ThisAlgo = Algo<GR, V, C>;
    using Index = unsigned int;
    using MCFArcMapV = typename GR::template ArcMap<V>;
    using MCFNodeMapV = typename GR::template NodeMap<V>;
    typedef typename ThisAlgo::ProblemType ProblemType;

    enum str_par_type_LEMON {
        strDMXFile = strLastParCDAS, ///< DMX filename to output the instance
        strLastParLEMON              ///< first allowed parameter value for derived classes
                                     /**< convenience value for easily allow derived classes
                                      * to further extend the set of types of return codes */
    };

    enum LEMON_sol_type { INFEASIBLE, OPTIMAL, UNBOUNDED };

    static constexpr int kErrorStatus = -1;

    MCFLemonSolver(void) {
        um = NULL;
        cm = NULL;
        bm = NULL;
    }

    ~MCFLemonSolver(void) {
        delete um;
        delete cm;
        delete bm;
    }

    /*--------------------------------------------------------------------------*/
    /*-------------------------- PUBLIC METHODS --------------------------------*/
    /*--------------------------------------------------------------------------*/
    /*-------------------------- OTHER INITIALIZATIONS -------------------------*/
    /*--------------------------------------------------------------------------*/
    /// register the Solver to a MCFBlock
    /** Provides the [MCF]Block encoding the problem instance that the
     * Solver will solve. It is entirely implemented in this class because this
     * is done uniformly for all LEMON algorithms. */

    void set_Block(Block *block) override;

    /*--------------------------------------------------------------------------*/
    /*--------------------- METHODS FOR SOLVING THE Block ----------------------*/
    /*--------------------------------------------------------------------------*/
    /// (try to) solve the MCF encoded in the MCFBlock
    /** Basically invokes the run() method of the underlying LEMON Algo. A lot
     * of the preparatory steps (locking the Block and the Solver, printing the
     * DMX file if required ...) are common to all the Algo and therefore are
     * implemented in this class; a guts_of_compute() method is inkoked at the
     * right time that eed be implemented in specialised classe. */

    int compute(bool changedvars = true) override;

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    /// @brief used for obtain default value of a string parameter
    /// @param par
    /// @return default valule of parameter par, if exists
    [[nodiscard]] const std::string &get_dflt_str_par(idx_type par) const override {
        if (par > strLastParLEMON) {
            throw std::invalid_argument("Invalid str parameter: out_of_range " + std::to_string(par));
        }
        static const std::string _empty;
        if (par == strLastParLEMON)
            return (_empty);

        return (CDASolver::get_dflt_str_par(par));
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    /// @brief used for get the value of string parameters
    /// @param par
    /// @return value of parameter indexed by par
    [[nodiscard]] const std::string &get_str_par(idx_type par) const override {
        if (par == strDMXFile)
            return (this->f_dmx_file);

        return (get_dflt_str_par(par));
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    /// @return number of str algorithimc parameters
    [[nodiscard]] idx_type get_num_str_par(void) const override { return (strLastParLEMON); }

    /*--------------------------------------------------------------------------*/
    /*---------------------- METHODS FOR READING RESULTS -----------------------*/
    /*--------------------------------------------------------------------------*/

    double get_elapsed_time(void) const override { return (this->ticks); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    OFValue get_lb(void) override { return OFValue(f_algo->totalCost()); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     * -*/

    OFValue get_ub(void) override { return OFValue(f_algo->totalCost()); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    int get_status(void) const { return (this->status); }

    /*--------------------------------------------------------------------------*/

    bool has_var_solution(void) override {
        switch (this->get_status()) {
        case (ThisAlgo::ProblemType::OPTIMAL):
        case (ThisAlgo::ProblemType::UNBOUNDED):
            return (true);
        default:
            return (false);
        }
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    bool has_dual_solution(void) override {
        switch (this->get_status()) {
        case (ThisAlgo::ProblemType::OPTIMAL):
        case (ThisAlgo::ProblemType::INFEASIBLE):
            return (true);
        default:
            return (false);
        }
    }

    void get_dual_solution(Configuration *solc = nullptr) override {
        auto MCFB = static_cast<MCFBlock *>(f_Block);
        Index i = 0;
        for (typename GR::NodeIt n(*dgp); n != INVALID; ++n) {
            MCFB->set_pi(i, f_algo->potential(n));
        }
    }

    /*--------------------------------------------------------------------------*/

    void get_var_solution(Configuration *solc = nullptr) override {
        auto MCFB = static_cast<MCFBlock *>(f_Block);
        Index i = 0;
        for (typename GR::ArcIt a(*dgp); a != INVALID; ++a) {
            MCFB->set_x(i, f_algo->flow(a));
        }
    }

    void get_var_direction(Configuration *dirc = nullptr) override {
        throw(std::logic_error("LEMONSolver:get_dual_direction() called"));
    }

    /*--------------------------------------------------------------------------*/
    /// returns false until we understand if and how LEMON does is

    bool has_var_direction(void) override { return (false); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     * -*/
    /// returns false until we understand if and how LEMON does is

    bool has_dual_direction(void) override { return (false); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     * -*/
    /// write the current dual direction in the Constraint of the MCFBlock
    /** Write the current unbounded dual direction, i.e., a cut separating two
     * shores to that the residual demand in one is greater than the capacity
     * across them, in the Constraint of the Block, in particular in the dual
     * variables of the flow conservation ones. To keep the same format as
     * MCFBlock::get_Solution() and MCFBlock::map[forward/back]_Modification(),
     * the Configuration *solc can be used to "partly" save it. In particular,
     * if solc != nullptr, it is a SimpleConfiguration< int >, and solc->f_value
     * == 1, then *nothing is done*, since the Configuration is meant to say
     * "only save/map the primal information". In all other cases, the direction
     * (cut) is saved.
     *
     * Or, rather, THIS SHOULD BE DONE, BUT THE METHOD IS NOT IMPLEMENTED yet.
     */

    void get_dual_direction(Configuration *dirc = nullptr) override {
        throw(std::logic_error("LEMONSolver:get_dual_direction() called"));
    }

    /*--------------------------------------------------------------------------*/
    /*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
    /*--------------------------------------------------------------------------*/

    /** @}
     * ---------------------------------------------------------------------*/
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

    void add_Modification(sp_Mod &mod) override;

    /*--------------------------------------------------------------------------*/
    /*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
    /*--------------------------------------------------------------------------*/
=======
/** @defgroup MCFLemonSolver_CLASSES Classes in MCFLemonSolver.h
 *  @{ */



  template< typename Algo >
  struct Fields {};
   
  enum LEMON_sol_type
  {
  UNSOLVED, //= NULL, ///< the problem has not been solved yet
  OPTIMAL,           ///< the problem has been solved
  KSTOPTIME, //= NULL,     ///< the problem has been stopped because of time limit
  INFEASIBLE,   ///< the problem is provably infeasible
  UNBOUNDED,    ///< the problem is provably unbounded
  KERROR //= NULL         ///< the problem has been stopped because of unrecoverable error
  };// end( sol_type )

/*--------------------------------------------------------------------------*/
/*--------------------------------- BASE CLASS -----------------------------*/
/*--------------------------------------------------------------------------*/
/** This class contains the fields and functions that every specialized version
 * of MCFLemonSolver<> has in common. For now this class contains:
 * 
 * set_Block() function
 * f_algo private field of type pointer to Algo<GR, V, C>
 * dgp private field of type pointer to GR
 * 
*/
template< template< typename , typename , typename > typename Algo ,
          typename GR , typename V , typename C >
  requires LEMONGraph< GR >
class MCFLemonSolverBase: virtual public CDASolver {

  public:

  enum str_par_type_LEMON_NS {
  strDMXFile = strLastParCDAS ,  ///< DMX filename to output the instance
  strLastParLEMON    ///< first allowed parameter value for derived classes
                   /**< convenience value for easily allow derived classes
                    * to further extend the set of types of return codes */
  };

  /// constructor: Initializes f_algo
  /** Void constructor. Initialize f_algo to nullptr
  */

  MCFLemonSolverBase(): CDASolver() {
    f_algo = NULL;
  }

  /// destructor: delete heap memory
  /**
   * Void destructor. delete f_algo and dgp fields from memory
   * Actually we have trouble with deleting dgp.
  */
  ~MCFLemonSolverBase( void ) {
    delete f_algo;
    delete dgp;
  
  }

  ///set_block function. Initialize dgp and f_algo
  /**
   * Prepare the block to be solved, initialize structure useful for LEMON algorithm
  */

  void set_Block(Block *block) override
    {
     
      
      if (block == f_Block) // actually doing nothing
        return;             // cowardly and silently return
      
      delete f_algo;
      f_algo=NULL;

      Solver::set_Block(block); // attach to the new Block

      if (block)
      { // this is not just resetting everything
        auto MCFB = dynamic_cast<MCFBlock *>(block);
        if (!MCFB)
          throw(std::invalid_argument(
              "MCFSolver:set_Block: block must be a MCFBlock"));

        bool owned = MCFB->is_owned_by(f_id);
        if ((!owned) && (!MCFB->read_lock()))
          throw(std::logic_error("cannot acquire read_lock on MCFBlock"));

        // create and clear new Graph (ListDigraph or SmartDigraph)
        dgp = new GR;
        dgp->clear();

        //Actually reserve get_MaxNNodes() node space in dgp
        dgp->reserveNode( MCFB->get_MaxNNodes() );
        MCFBlock::Index n = MCFB->get_NNodes();

        for( MCFBlock::Index i = 0; i < n; ++i)
          dgp->addNode();
        
        //Reserve gete_MaxNArcs() arcs space
        dgp->reserveArc( MCFB->get_MaxNArcs() );
        MCFBlock::Index m = MCFB->get_NArcs();

        //Get starting node subset and ending node subset
        MCFBlock::c_Subset & sn = MCFB->get_SN();
        MCFBlock::c_Subset & en = MCFB->get_EN();

        //Add arc in dgp
        /// sn[i] - 1 and en[i] - 1 are used because sn,en nodes start from 1
        for( MCFBlock::Index i = 0; i < m; ++i){
          dgp->addArc( dgp->nodeFromId( sn[i] - 1) , dgp->nodeFromId( en[i] - 1 ) );
        }


        //New instance of Lemon Algorithm
        f_algo = new Algo< GR , V, C >(*dgp);

        //Defining names for types for readability
        using MCFArcMapV = typename GR::template ArcMap< V >;
        using MCFNodeMapV = typename GR::template NodeMap< V >;


        //Now we are going to fill up ArcMap
        //This is the case of upperMap 
        if(!MCFB->get_U().empty())
        {
          MCFArcMapV um(*dgp);
          MCFBlock::c_Vec_FNumber & u = MCFB->get_U();
          for( MCFBlock::Index i = 0; i < m; ++i){
          um.set( dgp->arcFromId(i), u[i]);
          }
          f_algo->upperMap(um);
        }

        //This is the case of CostMap
        if(!MCFB->get_C().empty())
        {
          MCFArcMapV cm(*dgp);
          MCFBlock::c_Vec_FNumber & c = MCFB->get_C();
          for( MCFBlock::Index i = 0; i < m; ++i){
            cm.set( dgp->arcFromId(i), c[i]);
          }
          f_algo->costMap(cm);

        }

        //This is the case of supplyMap
        if(!MCFB->get_B().empty())
        {

          MCFNodeMapV bm(*dgp);
          MCFBlock::c_Vec_FNumber & b = MCFB->get_B();
          for( MCFBlock::Index i = 0; i < n; i++){
            bm.set( dgp->nodeFromId(i), -b[i]);
          }
          f_algo->supplyMap(bm);

        }

        
        // TODO: PreProcess() changes the internal data of the MCFSolver using
        //       information about how the data of the MCF is *now*. If the data
        //       changes, some of the deductions (say, reducing the capacity but
        //       of some arcs) may no longer be correct and they should be undone,
        //       there isn't any proper way to handle this. Thus, PreProcess() has
        //       to be disabled for now; maybe later on someone will take care to
        //       make this work (or maybe not).
        // MCFC::PreProcess();

        // once done, read_unlock the MCFBlock (if it was read-lock()-ed)
        if (!owned)
          MCFB->read_unlock();

        
      }
      
    } // end( set_Block )

>>>>>>> 3587a0f95544d63e767d768c535e345b48386907

  protected:

<<<<<<< HEAD
    void guts_of_constructor(void) { f_algo = nullptr; }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     * -*/

    void guts_of_destructor(void) {
        delete f_algo;
        dgp->clear();
        delete dgp;
    }

    /*--------------------------------------------------------------------------*/

    virtual void guts_of_compute(void) = 0;
=======
    Algo<GR, V, C> * f_algo; ///f_algo represents the algorithm used by Lemon for solving the MCFBlock
    GR *dgp; ///dgp represents the directed graph implemented by two classes by Lemon (ListDigraph and SmartDigraph)
>>>>>>> 3587a0f95544d63e767d768c535e345b48386907

}; // end( MCFLemonSolverBase )

<<<<<<< HEAD
    void guts_of_set_Block(MCFBlock *MCFB);

    /*--------------------------------------------------------------------------*/
    void process_outstanding_Modification(void);

    /*--------------------------------------------------------------------------*/
    /*---------------------------- PROTECTED FIELDS ---------------------------*/
    /*--------------------------------------------------------------------------*/

    std::string f_dmx_file;
    ///< string for DMX file output

    Algo<GR, V, C> *f_algo;
    ProblemType status = ProblemType::INFEASIBLE;

    GR *dgp;
    /**< represents the directed graph implemented by two classes by Lemon
     * (ListDigraph and SmartDigraph) */

    bool cost_changed = false;
    // represents when a cost Modification is executed, useful for updating
    // _cost vector in f_algo
    bool cap_changed = false;
    // represents when a capacity Modification is executed, useful for updating
    // _upper vector in f_algo
    bool supply_changed = false;
    // represents when a supply Modification is executed, useful for updating
    // _supply vector in f_algo

    int n_arcs;
    // number of arcs in the graph (closed arcs count like deleted)
    int n_arcs_added = 0;
    // number of arcs added in the graph with addArc/openArc
    int n_arcs_deleted = 0;
    // number of arcs deleted in the graph with rmvArc/closeArc

    /*--------------------------------------------------------------------------*/
    /*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
    /*--------------------------------------------------------------------------*/

  private:
    /*--------------------------------------------------------------------------*/
    /*--------------------------- PRIVATE FIELDS -------------------------------*/
    /*--------------------------------------------------------------------------*/

    double ticks; // Elapsed time in ticks for compute() method

    MCFArcMapV *um;
    // ArcMap that contains the upper bound on the flow of each arc

    MCFArcMapV *cm;
    // ArcMap that contains the cost of each arc

    MCFNodeMapV *bm;
    // NodeMap that contains the supply values of each node;

    /*--------------------------------------------------------------------------*/

}; // end( class MCFLemonSolver< Algo , GR , C , V > )

/*--------------------------------------------------------------------------*/
/*------------------------- SPECIALIZED CLASSES ----------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------- MCFLemonSolverNetworkSimplex ------------------------*/
/*--------------------------------------------------------------------------*/
/** Specialized MCFLemonSolver< NetworkSimplex , GR , V , C > that contains
 * specialized  compute() method, enums for indexing algorithimc parameters
 * and function  set/get_*_par for manage them.
 *
 *  Template parameters are:
 *
 *  - NetworkSimplex, implements the primal Network Simplex algorithm for
 *    finding a minimum cost flow. This algorithm is a highly efficient
 *    specialized version of the linear programming simplex method directly
 *    for the minimum cost flow problem.
 *
 *  - GR that represents the directed graph, the possibilities are described
 *    in the file MCFLemonSolver.h
 *
 *  - V, which is the type of flows / deficits; typically, double can be used
 *    for maximum compatibility, but long (or even smaller) would yeld better
 *    performances;
 *
 * - C, which is the type of ar costs; typically, double can be used for
 *    maximum compatibility, but long (or even smaller) would yeld better
 *    performances; */

template <typename GR, typename V, typename C>
class MCFLemonSolverNetworkSimplex : public MCFLemonSolver<NetworkSimplex, GR, V, C> {
    /*--------------------------------------------------------------------------*/
    /*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
    /*--------------------------------------------------------------------------*/

  public:
    /*--------------------------------------------------------------------------*/
    /*---------------------------- PUBLIC TYPES --------------------------------*/
    /*--------------------------------------------------------------------------*/
    // large batch of "using", unfortunately needed due to the fact that on
    // the first pass of compiling a template, where only non-dependent names
    // are looked up, the compiler does not "see" the base class MCFLemonSolver
    // and all its ancestors (CDASolver, Solver, ThinComputeInterface) and all
    // their names

    using BaseClass = MCFLemonSolver<NetworkSimplex, GR, V, C>;

    using BaseClass::f_algo;
    using BaseClass::intLastParCDAS;
    using BaseClass::status;
    using BaseClass::strLastParLEMON;
    using BaseClass::str_par_type_LEMON::strDMXFile;
    using SMSpp_di_unipi_it::CDASolver::kBlockLocked;
    using SMSpp_di_unipi_it::CDASolver::kInfeasible;
    using SMSpp_di_unipi_it::CDASolver::kStopTime;
    using SMSpp_di_unipi_it::Solver::f_Block;
    using SMSpp_di_unipi_it::Solver::f_id;
    using SMSpp_di_unipi_it::Solver::lock;
    using SMSpp_di_unipi_it::Solver::OFValue;
    using SMSpp_di_unipi_it::Solver::unlock;
    using SMSpp_di_unipi_it::ThinComputeInterface::idx_type;
    using SMSpp_di_unipi_it::ThinComputeInterface::kUnEval;

    using typename BaseClass::ThisAlgo;
    using NSPivotRule = typename NetworkSimplex<GR, V, C>::PivotRule;

    /*--------------------------------------------------------------------------*/
    // enums for handling the extra parameters

    enum LEMON_NS_dbl_par_type {
        dblLastParLEMON_NS ///< first allowed parameter value for derived
                           ///< classes
        /**< convenience value for easily allow derived classes
         * to further extend the set of types of return codes */
    };

    enum LEMON_NS_int_par_type {
        kPivot = intLastParCDAS,
        intLastParLEMON_NS ///< first allowed parameter value for derived
                           ///< classes
        /**< convenience value for easily allow derived classes
         * to further extend the set of types of return codes */
    };

    /** @}
     * ---------------------------------------------------------------------*/
    /*------- CONSTRUCTING AND DESTRUCTING MCFLemonSolverNetworkSimplex --------*/
    /*--------------------------------------------------------------------------*/
    /** @name Constructing and destructing MCFLemonSolverNetworkSimplex
     *  @{ */

    /// constructor: Initializes algorithm parameters
    /** Void constructor. Define f_pivot_rule to the default algorithm
     * parameters used by NetworkSimplex.
     */

    MCFLemonSolverNetworkSimplex(void) {
        BaseClass::guts_of_constructor();
        f_pivot_rule = NSPivotRule::BLOCK_SEARCH;
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    /// destructor: delete algorithm parameter
    /** Does nothing special, delete Fields f_pivot_rule, an algorithmic
     * parameter of  NetworkSimplex */

    ~MCFLemonSolverNetworkSimplex() { BaseClass::guts_of_destructor(); }

    /*--------------------------------------------------------------------------*/
    /// @brief set the parameter par with value
    /// @param par
    /// @param value

    void set_par(idx_type par, int value) override {
        if (par == kPivot) {
            if ((value < 0) || (value > 4))
                throw(std::invalid_argument("Error: invalid kPivot " + std::to_string(value)));

            if (value == f_pivot_rule)
                return; // nothing is changed

            f_pivot_rule = NSPivotRule(value);
            return;
        }

        CDASolver::set_par(par, value);
        return;
    }

    /** @}
     * ---------------------------------------------------------------------*/
    /*--------------------- METHODS FOR SOLVING THE Block ----------------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Solving the MCF encoded by the current MCFBlock
     *  @{ */

    /// (try to) solve the MCF encoded in the MCFBlock

    /*--------------------------------------------------------------------------*/
    /*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
    /*--------------------------------------------------------------------------*/

    ///@return number of int algorithmic parameters
    [[nodiscard]] idx_type get_num_int_par(void) const override { return (intLastParLEMON_NS); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     * -*/

    /// @return number of dbl algorithmic parameters
    [[nodiscard]] idx_type get_num_dbl_par(void) const override { return (dblLastParLEMON_NS); }

    /*--------------------------------------------------------------------------*/
    /// @brief used for obtain default value of an int parameter
    /// @param par
    /// @return default value of parameter par, if exists
    [[nodiscard]] int get_dflt_int_par(idx_type par) const override {
        if (par > intLastParLEMON_NS)
            throw std::invalid_argument("Invalid int parameter: out_of_range " + std::to_string(par));

        if (par == kPivot)
            return (NSPivotRule::BLOCK_SEARCH);

        return (CDASolver::get_dflt_int_par(par));
    }

    /*--------------------------------------------------------------------------*/
    /// @brief used for get the value of int parameters
    /// @param par
    /// @return value of parameter indexed by par

    [[nodiscard]] int get_int_par(idx_type par) const override {
        if (par == kPivot)
            return (f_pivot_rule);

        return (get_dflt_int_par(par));
    }

    /*--------------------------------------------------------------------------*/
    /// @brief used for convert string to int parameter's index
    /// @param name
    /// @return an index that denotes parameter name

    [[nodiscard]] idx_type int_par_str2idx(const std::string &name) const override {
        if (name == "kPivot")
            return (kPivot);

        return (CDASolver::int_par_str2idx(name));
    }

    /*--------------------------------------------------------------------------*/
    /// @brief used for convert int index idx to phrasal rapresentation
    /// @param idx
    /// @return a string that denotes index idx parameter

    [[nodiscard]] const std::string &int_par_idx2str(idx_type idx) const override {
        if (idx > intLastParLEMON_NS)
            throw(std::invalid_argument("int_par_idx2str: index out_of_range " + std::to_string(idx)));

        static const std::string par = "kPivot";
        if (idx == kPivot)
            return (par);

        return (CDASolver::int_par_idx2str(idx));
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
    /// @brief used for convert dbl index idx to phrasal rapresentation
    /// @param idx
    /// @return a string that denotes index idx parameter

    [[nodiscard]] const std::string &dbl_par_idx2str(idx_type idx) const override {
        if (idx > dblLastParLEMON_NS) {
            throw std::invalid_argument("Invalid index: out_of_range " + std::to_string(idx));
        }

        switch (idx) {
        default:
            break;
        }

        return (CDASolver::dbl_par_idx2str(idx));
    }

    void guts_of_compute(void) override { status = f_algo->run(NSPivotRule(f_pivot_rule)); }

    /** @}
     * ---------------------------------------------------------------------*/
    /*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
    /*--------------------------------------------------------------------------*/

  protected:
    /*--------------------------------------------------------------------------*/
    /*---------------------------- PROTECTED FIELDS ---------------------------*/
    /*--------------------------------------------------------------------------*/

    NSPivotRule f_pivot_rule;

    /*--------------------------------------------------------------------------*/
    /*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
    /*--------------------------------------------------------------------------*/

  private:
    /*--------------------------------------------------------------------------*/
    /*-------------------------- PRIVATE METHODS -------------------------------*/
    /*--------------------------------------------------------------------------*/

    SMSpp_insert_in_factory_h;

    /*--------------------------------------------------------------------------*/

}; // end( class MCFLemonSolver<NetworkSimplex, GR, V, C>  Specialization)

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------------- CYCLECANCELING --------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/** Specialized MCFLemonSolverCycleCanceling< GR, V, C> that contains
 *  guts_of_compute() method, enums for indexing algorithimc parameters and
 * function set/get_*_par for manage them.
 *
 *  Template parameters are:
 *
 *  -  CycleCanceling implements three different cycle-canceling algorithms for
 * finding a minimum cost flow. The most efficent one is the Cancel-and-tighten
 * algorithm, thus it is the default method. It runs in strongly polynomial
 * time, but in practice, it is typically orders of magnitude slower than the
 * scaling algorithms and NetworkSimplex.
 *
 *  - GR that represents the directed graph, the possibilities are described in
 * the file MCFLemonSolver.h at line 339
 *
 *  - V, which is the type of flows / deficits; typically, double can be used
 *    for maximum compatibility, but long (or even smaller) would yeld better
 *    performances;
 *
 * - C, which is the type of ar costs; typically, double can be used for
 *    maximum compatibility, but long (or even smaller) would yeld better
 *    performances;
 */
template <typename GR, typename V, typename C>
class MCFLemonSolverCycleCanceling : public MCFLemonSolver<CycleCanceling, GR, V, C> {
    /*--------------------------------------------------------------------------*/
=======
/*--------------------------------------------------------------------------*/
/*------------------------ CLASS MCFLemonSolver ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// CDASolver for MCFBlock based on the LEMON project
/** The MCFLemonSolver implements Solver interface for MCFBlock that represent
 * (Linear) Min-Cost Flow (MCF) problems, using algorithms of LEMON library.
 * Because MCF is a Linear Program it has a(n exact) dual, and therefore
 * MCFLemonSolver implements the CDASolver interface for also giving out dual
 * information.
 *
 * The MCFLemonSolver is template over four different types:
 *
 * - GR, which is the type of graph:
 * 
 *   The possibilities are:
 *   - SmartDigraph is a simple and fast directed graph implementation
 *     It is quite memory efficient but at the price that it does not support
 *     node and arc deletion.
 *     It will be hard find the best way for modification.
 * 
 *   - ListDigraph, a versatile and fast directed graph implementation based
 *     on linked lists that are stored in std::vector structures.
 *     This class provides only linear time counting for nodes and arcs.
 *     It support node and arc deletion, useful for modification.
 *
 * - V, which is the type of flows / deficits; typically, double can be used
 *   for maximum compatibility, but int (or even smaller) would yeld better
 *   performances;
 *
 * - C, which is the type of ar costs; typically, double can be used for
 *   maximum compatibility, but int (or even smaller) would yeld better
 *   performances;
 *
 * - Algo, which is the specific algorithm (itself, template over GR, V, and
 *   C) implemented in the LEMON class:
 *   The possibilities are:
 * 
 *   - NetworkSimplex implements the primal Network Simplex algorithm for finding
 *     a minimum cost flow. This algorithm is a highly efficient specialized version
 *     of the linear programming simplex method directly for the minimum cost flow
 *     problem.
 *     
 *   - CycleCanceling implements three different cycle-canceling algorithms for finding 
 *     a minimum cost flow. The most efficent one is the Cancel-and-tighten algorithm,
 *     thus it is the default method. It runs in strongly polynomial time, but in practice, 
 *     it is typically orders of magnitude slower than the scaling algorithms and NetworkSimplex.
 * 
 *   - CostScaling implements a cost scaling algorithm that performs push/augment and
 *     relabel operations for finding a minimum cost flow. It is a highly efficient primal-dual
 *     solution method, which can be viewed as the generalization of the preflow push-relabel
 *     algorithm for the maximum flow problem. It is a polynomial algorithm.
 * 
 *   - CapacityScaling implements the capacity scaling version of the successive shortest path
 *     algorithm for finding a minimum cost flow. It is an efficient dual solution method,
 *     which runs in polynomial time.
 *     In special case it can be more efficient than CostScaling and NetworkSimplex algorithms.
 * 
 *   In general, NetworkSimplex and CostScaling are the fastest implementations available in LEMON
 *   for solving this problem.
 *    
 *   Note that scaling-type algorithms may behave in different ways according
 *   to which combination of V and C is used, and there are different
 *   "traits" for this which are another scaling parameter. However, in order
 *   to make MCFLemonSolver class template over always the same number of
 *   template parameters we fix the use of the default trait, which is why
 *   SMSppCapacityScaling and SMSppCostScaling are defined (as template over
 *   < GR , V , C >) that are meant to be used instead of the original
 *   CapacityScaling and CostScaling. */

template< template< typename , typename , typename > typename Algo ,
          typename GR , typename V , typename C >
  requires LEMONGraph< GR >
class MCFLemonSolver: virtual public CDASolver, Fields< Algo<GR, V, C> > {
};
          

  /*--------------------------------------------------------------------------*/
  /*------------------------- CLASS MCFLemonState ---------------------------*/
  /*--------------------------------------------------------------------------*/
  /// class to describe the "internal state" of a MCFLemonSolver
  /** Derived class from State to describe the "internal state" of a MCFLemonSolver,
   *  i.e., a MCFClass::MCFState (*). Since MCFClass::MCFState does not allow
   *  serialization, all that part does not work.  */

  class MCFLemonState : public State
  {
>>>>>>> 3587a0f95544d63e767d768c535e345b48386907
    /*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
    /*--------------------------------------------------------------------------*/

    using BaseClass = MCFLemonSolver<CycleCanceling, GR, V, C>;

    using BaseClass::intLastParCDAS;
    // using BaseClass::idx_type;

    using BaseClass::dgp;
    using BaseClass::f_algo;
    using BaseClass::status;
    using BaseClass::strLastParLEMON;
    using BaseClass::str_par_type_LEMON::strDMXFile;
    using SMSpp_di_unipi_it::CDASolver::dblLastParCDAS;
    using SMSpp_di_unipi_it::CDASolver::kBlockLocked;
    using SMSpp_di_unipi_it::CDASolver::kInfeasible;
    using SMSpp_di_unipi_it::CDASolver::kStopTime;
    using SMSpp_di_unipi_it::Solver::f_Block;
    using SMSpp_di_unipi_it::Solver::f_id;
    using SMSpp_di_unipi_it::Solver::lock;
    using SMSpp_di_unipi_it::Solver::OFValue;
    using SMSpp_di_unipi_it::Solver::unlock;
    using SMSpp_di_unipi_it::ThinComputeInterface::idx_type;
    using SMSpp_di_unipi_it::ThinComputeInterface::kUnEval;
    using typename BaseClass::ThisAlgo;
    using CCMethod = typename ThisAlgo::Method;

  public:
<<<<<<< HEAD
    /*--------------------------------------------------------------------------*/
    /*---------------------------- PUBLIC TYPES --------------------------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Public Types
     *  @{ */
    const int kErrorStatus = -1;
    /*
     * ---------------------------------------------------------------------*/
    /*-------------- CONSTRUCTING AND DESTRUCTING MCFLemonSolver
     * ---------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Constructing and destructing MCFLemonSolver
     *  @{ */

    /// constructor: assign the default parameter
    /** Void constructor: Build f_method with default parameter */

    MCFLemonSolverCycleCanceling(void) : MCFLemonSolver<CycleCanceling, GR, V, C>() {
        BaseClass::guts_of_constructor();
        f_method = CycleCanceling<GR, V, C>::Method::CANCEL_AND_TIGHTEN;
    }

    /// destructor:
    /**Void destructor: Delete f_method from memory */
    ~MCFLemonSolverCycleCanceling(void) { BaseClass::guts_of_destructor(); }

    /*--------------------------------------------------------------------------*/

    enum LEMON_CC_int_par_type {
        kMethod = intLastParCDAS,
        intLastParLEMON_CC ///< first allowed parameter value for derived
                           ///< classes
        /**< convenience value for easily allow derived classes
         * to further extend the set of types of return codes */
    };

    /*--------------------------------------------------------------------------*/

    /*--------------------------------------------------------------------------*/

    /// public enum for the type of the solution
    /// to MCFLemonSolver< CycleCanceling , C , V >

    enum dbl_par_type_LEMON_CC {
        dblCycleCancelingFactor = dblLastParCDAS,
        dblLastParLEMON_CC ///< first allowed parameter value for derived
                           ///< classes
        /**< convenience value for easily allow derived classes
         * to further extend the set of types of return codes */
    };

    /*--------------------------------------------------------------------------*/
    /// @brief set the parameter par with value
    /// @param par
    /// @param value
    void set_par(idx_type par, int value) override {

        if (par == kMethod) {
            if ((value < 0) || (value > 4))
                throw(std::invalid_argument("Error: invalid kMethod " + std::to_string(value)));

            if (value == f_method)
                return; // nothing is changed

            f_method = CCMethod(value);
            return;
        }

        CDASolver::set_par(par, value);
        return;
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    /*--------------------------------------------------------------------------*/
    /*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
    /*--------------------------------------------------------------------------*/

    ///@return number of int algorithmic parameters
    [[nodiscard]] idx_type get_num_int_par(void) const override { return (intLastParLEMON_CC); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    /// @return number of dbl algorithmic parameters
    [[nodiscard]] idx_type get_num_dbl_par(void) const override { return (dblLastParLEMON_CC); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    /*--------------------------------------------------------------------------*/

    /// @brief used for obtain default value of an int parameter
    /// @param par
    /// @return default value of parameter par, if exists
    [[nodiscard]] int get_dflt_int_par(idx_type par) const override {
        if (par > intLastParLEMON_CC) {
            throw(std::invalid_argument(std::to_string(par)));
        }

        switch (par) {
        case kMethod:
            return CycleCanceling<GR, int, int>::Method::CANCEL_AND_TIGHTEN;
        default:
            return (CDASolver::get_dflt_int_par(par));
        }
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    /// @brief used for obtain default value of a double parameter
    /// @param par
    /// @return default value of parameter par, if exists
    [[nodiscard]] double get_dflt_dbl_par(idx_type par) const override {
        if (par > dblLastParLEMON_CC) {
            throw(std::invalid_argument(std::to_string(par)));
        }

        switch (par) {
        default:
            return (CDASolver::get_dflt_dbl_par(par));
        }
    }

    /*--------------------------------------------------------------------------*/

    /// @brief used for get the value of int parameters
    /// @param par
    /// @return value of parameter indexed by par
    [[nodiscard]] int get_int_par(idx_type par) const override {

        if (par == kMethod) {
            return f_method;
        }

        return (get_dflt_int_par(par));
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    /// @brief used for get the value of dbl parameters
    /// @param par
    /// @return value of parameter indexed by par
    [[nodiscard]] double get_dbl_par(idx_type par) const override {
        // Da finire parametri algoritmici dbl
        return (get_dflt_dbl_par(par));
    }

    /*--------------------------------------------------------------------------*/

    /// @brief used for convert string to int parameter's index
    /// @param name
    /// @return an index that denotes parameter name
    [[nodiscard]] idx_type int_par_str2idx(const std::string &name) const override {
        if (name == "kMethod")
            return (kMethod);

        return (CDASolver::int_par_str2idx(name));
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    /// @brief used for convert string to dbl parameter's index
    /// @param name
    /// @return an index that denotes parameter name
    [[nodiscard]] idx_type dbl_par_str2idx(const std::string &name) const override {
        return (CDASolver::dbl_par_str2idx(name));
    }

    /*--------------------------------------------------------------------------*/

    /// @brief used for convert int index idx to phrasal rapresentation
    /// @param idx
    /// @return a string that denotes index idx parameter
    [[nodiscard]] const std::string &int_par_idx2str(idx_type idx) const override {

        if (idx > intLastParLEMON_CC) {
            throw std::invalid_argument(std::to_string(idx));
        }
        static const std::string par = "kMethod";
        switch (idx) {
        case kMethod:
            return (par);
        default:
            break;
        }

        return (CDASolver::int_par_idx2str(idx));
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    /// @brief used for convert dbl index idx to phrasal rapresentation
    /// @param idx
    /// @return a string that denotes index idx parameter
    [[nodiscard]] const std::string &dbl_par_idx2str(idx_type idx) const override {
        if (idx > dblLastParLEMON_CC) {
            throw std::invalid_argument(std::to_string(idx));
        }

        switch (idx) {
        default:
            break;
        }

        return (CDASolver::dbl_par_idx2str(idx));
    }

    /*--------------------------------------------------------------------------*/
    /*-------------------------------- FRIENDS ---------------------------------*/
    /*--------------------------------------------------------------------------*/

    friend class MCFLemonState; // make MCFSolverState friend

    /** @}
     * ---------------------------------------------------------------------*/
    /*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
    /*--------------------------------------------------------------------------*/

  protected:
    /*--------------------------------------------------------------------------*/
    /*-------------------------- PROTECTED METHODS -----------------------------*/
    /*--------------------------------------------------------------------------*/

    void guts_of_compute() override { status = f_algo->run(CCMethod(f_method)); }

    /*--------------------------------------------------------------------------*/
    /*---------------------------- PROTECTED FIELDS ---------------------------*/
    /*--------------------------------------------------------------------------*/

    std::string f_dmx_file;
    std::string f_LMN_file = "lmn.txt";
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
    CCMethod f_method;

    double ticks; // Elaped time in ticks for compute() method
    /*--------------------------------------------------------------------------*/

}; // end( class MCFLemonSolver<CycleCanceling> specialization )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------------- CAPACITYSCALING -------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/** Specialized MCFLemonSolverCapacityScaling that contains guts_of_compute()
 * enums for indexing algorithimc parameters and function set/get_*_par for
 * manage them.
 *
 *  Template parameters are:
 *
 *   - CapacityScaling implements the capacity scaling version of the successive
 * shortest path algorithm for finding a minimum cost flow. It is an efficient
 * dual solution method, which runs in polynomial time. In special case it can
 * be more efficient than CostScaling and NetworkSimplex algorithms.
 *
 *  - GR that represents the directed graph, the possibilities are described in
 * the file MCFLemonSolver.h at line 339
 *
 *  - V, which is the type of flows / deficits; typically, l can be used
 *    for maximum compatibility, but long (or even smaller) would yeld better
 *    performances;
 *
 * - C, which is the type of ar costs; typically, double can be used for
 *    maximum compatibility, but long (or even smaller) would yeld better
 *    performances;
 */
template <typename GR, typename V, typename C>
class MCFLemonSolverCapacityScaling : public MCFLemonSolver<SMSppCapacityScaling, GR, V, C> {
    /*--------------------------------------------------------------------------*/
    /*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
    /*--------------------------------------------------------------------------*/

    using BaseClass = MCFLemonSolver<SMSppCapacityScaling, GR, V, C>;

    using BaseClass::intLastParCDAS;
    // using BaseClass::idx_type;

    using BaseClass::dgp;
    using BaseClass::f_algo;
    using BaseClass::status;
    using BaseClass::strLastParLEMON;
    using BaseClass::str_par_type_LEMON::strDMXFile;
    using SMSpp_di_unipi_it::CDASolver::kBlockLocked;
    using SMSpp_di_unipi_it::CDASolver::kInfeasible;
    using SMSpp_di_unipi_it::CDASolver::kStopTime;
    using SMSpp_di_unipi_it::Solver::f_Block;
    using SMSpp_di_unipi_it::Solver::f_id;
    using SMSpp_di_unipi_it::Solver::lock;
    using SMSpp_di_unipi_it::Solver::OFValue;
    using SMSpp_di_unipi_it::Solver::unlock;
    using SMSpp_di_unipi_it::ThinComputeInterface::idx_type;
    using SMSpp_di_unipi_it::ThinComputeInterface::kUnEval;

  public:
    /*--------------------------------------------------------------------------*/
    /*---------------------------- PUBLIC TYPES -----------------------------*/
    /*--------------------------------------------------------------------------*/
    const int kErrorStatus = -1;
    /*  *
     * ---------------------------------------------------------------------*/
    /*-------------- CONSTRUCTING AND DESTRUCTING MCFLemonSolver --------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Constructing and destructing MCFLemonSolver
     *  @{ */

    /// constructor: does nothing
    MCFLemonSolverCapacityScaling(void) : BaseClass() { BaseClass::guts_of_constructor(); }

    /// destructor: doesnothing
    ~MCFLemonSolverCapacityScaling(void) { BaseClass::guts_of_destructor(); }

    /*--------------------------------------------------------------------------*/

    enum LEMON_CS_int_par_type {
        intLastParLEMON_CS //< first allowed parameter value for derived classes
        /**< convenience value for easily allow derived classes
         * to further extend the set of types of return codes */
    };

    /*--------------------------------------------------------------------------*/

    enum LEMON_CS_dbl_par_type {
        dblLastParLEMON_CS //< first allowed parameter value for derived classes
        /**< convenience value for easily allow derived classes
         * to further extend the set of types of return codes */
    };

    /*--------------------------------------------------------------------------*/

    void set_par(idx_type par, int value) override {

        CDASolver::set_par(par, value);
        return;
    }

    /*--------------------------------------------------------------------------*/

    /*--------------------------------------------------------------------------*/
    /*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Handling the parameters of the MCFLemonSolver
     *
     * Each MCFLemonSolver< Algo > may have its own extra int / double
     * parameters. If this is the case, it will have to specialize the following
     * methods to handle them. The general definition just handles the case of
     * the
     *
     * intLastParCDAS ==> kReopt             whether or not to reoptimize
     *
     *  @{ */

    ///@return number of int algorithmic parameters
    [[nodiscard]] idx_type get_num_int_par(void) const override { return (intLastParLEMON_CS); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  -*/

    /// @return number of dbl algorithmic parameters
    [[nodiscard]] idx_type get_num_dbl_par(void) const override { return (dblLastParLEMON_CS); }

    /*--------------------------------------------------------------------------*/

    /// @brief used for obtain default value of an int parameter
    /// @param par
    /// @return default value of parameter par, if exists
    [[nodiscard]] int get_dflt_int_par(idx_type par) const override {
        if (par > intLastParLEMON_CS) {
            throw std::invalid_argument(std::to_string(par));
        }

        switch (par) {
        default:
            return (CDASolver::get_dflt_int_par(par));
        }
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     * -*/

    /// @brief used for obtain default value of an dbl parameter
    /// @param par
    /// @return default value of parameter par, if exists
    [[nodiscard]] double get_dflt_dbl_par(idx_type par) const override {
        if (par > dblLastParLEMON_CS) {
            throw std::invalid_argument(std::to_string(par));
        }

        switch (par) {
        default:
            return (CDASolver::get_dflt_dbl_par(par));
        }
    }

    /*--------------------------------------------------------------------------*/

    /// @brief used for get the value of int parameters
    /// @param par
    /// @return value of parameter indexed by par
    [[nodiscard]] int get_int_par(idx_type par) const override {
        // No int parameters
        return (get_dflt_int_par(par));
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    /// @brief used for get the value of dbl parameters
    /// @param par
    /// @return value of parameter indexed by par
    [[nodiscard]] double get_dbl_par(idx_type par) const override {
        // No dbl parameters
        return (get_dflt_dbl_par(par));
    }

    /*--------------------------------------------------------------------------*/

    /// @brief used for convert string to int parameter's index
    /// @param name
    /// @return an index that denotes parameter name
    [[nodiscard]] idx_type int_par_str2idx(const std::string &name) const override {
        return (CDASolver::int_par_str2idx(name));
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    /// @brief used for convert string to dbl parameter's index
    /// @param name
    /// @return an index that denotes parameter name
    [[nodiscard]] idx_type dbl_par_str2idx(const std::string &name) const override {
        return (CDASolver::dbl_par_str2idx(name));
    }

    /*--------------------------------------------------------------------------*/

    /// @brief used for convert int index idx to phrasal rapresentation
    /// @param idx
    /// @return a string that denotes index idx parameter
    [[nodiscard]] const std::string &int_par_idx2str(idx_type idx) const override {
        if (idx > intLastParLEMON_CS) {
            throw std::invalid_argument(std::to_string(idx));
        }

        switch (idx) {
        default:
            break;
        }

        return (CDASolver::int_par_idx2str(idx));
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    /// @brief used for convert dbl index idx to phrasal rapresentation
    /// @param idx
    /// @return a string that denotes index idx parameter
    [[nodiscard]] const std::string &dbl_par_idx2str(idx_type idx) const override {
        if (idx > dblLastParLEMON_CS) {
            throw std::invalid_argument(std::to_string(idx));
        }

        switch (idx) {
        default:
            break;
        }

        return (CDASolver::dbl_par_idx2str(idx));
    }

    /*--------------------------------------------------------------------------*/
    /*-------------------------------- FRIENDS ---------------------------------*/
    /*--------------------------------------------------------------------------*/

    friend class MCFLemonState; // make MCFSolverState friend

    /** @}
     * ---------------------------------------------------------------------*/
    /*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
    /*--------------------------------------------------------------------------*/

  protected:
    /*--------------------------------------------------------------------------*/
    /*-------------------------- PROTECTED METHODS -----------------------------*/
    /*--------------------------------------------------------------------------*/

    void guts_of_compute() override { status = f_algo->run(); }

    /*--------------------------------------------------------------------------*/
    /*---------------------------- PROTECTED FIELDS ---------------------------*/
    /*--------------------------------------------------------------------------*/

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

    double ticks; // Elaped time in ticks for compute() method
    /*--------------------------------------------------------------------------*/

}; // end( class MCFLemonSolver<SMSppCapacityScaling, GR, V, C>  Specialization)

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------------- COSTSCALING -----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/** Specialized MCFLemonSolverCostScaling< GR, V, C> that contains specialized
 *  compute() method, enums for indexing algorithimc parameters and function
 *  set/get_*_par for manage them.
 *
 *  Template parameters are:
 *
 *  - GR that represents the directed graph
 *
 *  - V, which is the type of flows / deficits; typically, double can be used
 *    for maximum compatibility, but long (or even smaller) would yeld better
 *    performances;
 *
 * - C, which is the type of ar costs; typically, double can be used for
 *    maximum compatibility, but long (or even smaller) would yeld better
 *    performances;
 */

template <typename GR, typename V, typename C>
class MCFLemonSolverCostScaling : public MCFLemonSolver<SMSppCostScaling, GR, V, C> {
    /*--------------------------------------------------------------------------*/
    /*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
    /*--------------------------------------------------------------------------*/

  public:
    using BaseClass = MCFLemonSolver<SMSppCostScaling, GR, V, C>;

    using BaseClass::intLastParCDAS;
    // using BaseClass::idx_type;

    using BaseClass::dgp;
    using BaseClass::f_algo;
    using BaseClass::status;
    using BaseClass::strLastParLEMON;
    using BaseClass::str_par_type_LEMON::strDMXFile;
    using SMSpp_di_unipi_it::CDASolver::kBlockLocked;
    using SMSpp_di_unipi_it::CDASolver::kInfeasible;
    using SMSpp_di_unipi_it::CDASolver::kStopTime;
    using SMSpp_di_unipi_it::Solver::f_Block;
    using SMSpp_di_unipi_it::Solver::f_id;
    using SMSpp_di_unipi_it::Solver::lock;
    using SMSpp_di_unipi_it::Solver::OFValue;
    using SMSpp_di_unipi_it::Solver::unlock;
    using SMSpp_di_unipi_it::ThinComputeInterface::idx_type;
    using SMSpp_di_unipi_it::ThinComputeInterface::kUnEval;

    using typename BaseClass::ThisAlgo;
    using CSMethod = typename ThisAlgo::Method;

    /*--------------------------------------------------------------------------*/
    /*---------------------------- PUBLIC TYPES --------------------------------*/
    /*--------------------------------------------------------------------------*/
    const int kErrorStatus = -1;

    /** @}
     * ---------------------------------------------------------------------*/
    /*-------------- CONSTRUCTING AND DESTRUCTING MCFLemonSolver ---------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Constructing and destructing MCFLemonSolver
     *  @{ */

    /// constructor: Build f_method
    /** Void constructor: Initialize f_method to default algorithmic parameter
     */

    MCFLemonSolverCostScaling(void) {
        BaseClass::guts_of_constructor();
        f_method = SMSppCostScaling<GR, V, C>::Method::PARTIAL_AUGMENT;
    }

    /// @brief delete f_method algorithmic parameter pointer
    /// @param
    ~MCFLemonSolverCostScaling(void) { BaseClass::guts_of_constructor(); }

    /*--------------------------------------------------------------------------*/

    enum LEMON_CS_int_par_type {
        kMethod = intLastParCDAS,
        intLastParLEMON_CS //< first allowed parameter value for derived classes
        /**< convenience value for easily allow derived classes
         * to further extend the set of types of return codes */
    };

    /*--------------------------------------------------------------------------*/
    enum LEMON_CS_dbl_par_type {
        dblLastParLEMON_CS //< first allowed parameter value for derived classes
        /**< convenience value for easily allow derived classes
         * to further extend the set of types of return codes */
    };

    /// @brief set the parameter par with value
    /// @param par
    /// @param value
    void set_par(idx_type par, int value) override {

        if (par == kMethod) {
            if ((value < 0) || (value > 4))
                throw(std::invalid_argument("Error: invalid kPivot " + std::to_string(value)));

            if (value == f_method)
                return; // nothing is changed

            f_method = CSMethod(value);
            return;
        }

        CDASolver::set_par(par, value);
        return;
    }

    /*--------------------------------------------------------------------------*/
    /*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
    /*--------------------------------------------------------------------------*/
    /** @name Handling the parameters of the MCFLemonSolver
     *
     * Each MCFLemonSolver< Algo > may have its own extra int / double
     * parameters. If this is the case, it will have to specialize the following
     * methods to handle them. The general definition just handles the case of
     * the
     *
     * intLastParCDAS ==> kReopt             whether or not to reoptimize
     *
     *  @{ */

    ///@return number of int algorithmic parameters
    [[nodiscard]] idx_type get_num_int_par(void) const override { return (intLastParLEMON_CS); }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    ///@return number of dbl algorithmic parameters
    [[nodiscard]] idx_type get_num_dbl_par(void) const override { return (dblLastParLEMON_CS); }

    /*--------------------------------------------------------------------------*/

    /// @brief used for obtain default value of an int parameter
    /// @param par
    /// @return default value of parameter par, if exists
    [[nodiscard]] int get_dflt_int_par(idx_type par) const override {

        if (par > intLastParLEMON_CS) {
            throw std::invalid_argument(std::to_string(par));
        }

        switch (par) {
        case kMethod:
            return SMSppCostScaling<GR, V, C>::Method::PARTIAL_AUGMENT;
        default:
            return (CDASolver::get_dflt_int_par(par));
        }
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    /// @brief used for obtain default value of an dbl parameter
    /// @param par
    /// @return default value of parameter par, if exists
    [[nodiscard]] double get_dflt_dbl_par(idx_type par) const override {
        if (par > intLastParLEMON_CS) {
            throw std::invalid_argument(std::to_string(par));
        }

        switch (par) {
        default:
            return (CDASolver::get_dflt_dbl_par(par));
        }
    }

    /*--------------------------------------------------------------------------*/

    /// @brief used for get the value of int parameters
    /// @param par
    /// @return value of parameter indexed by par
    [[nodiscard]] int get_int_par(idx_type par) const override {

        if (par == kMethod) {
            return f_method;
        }
        return (get_dflt_int_par(par));
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    /// @brief used for get the value of dbl parameters
    /// @param par
    /// @return value of parameter indexed by par
    [[nodiscard]] double get_dbl_par(idx_type par) const override { return (get_dflt_dbl_par(par)); }

    /*--------------------------------------------------------------------------*/

    /// @brief used for convert string to int parameter's index
    /// @param name
    /// @return an index that denotes parameter name
    [[nodiscard]] idx_type int_par_str2idx(const std::string &name) const override {
        if (name == "kMethod")
            return (kMethod);

        return (CDASolver::int_par_str2idx(name));
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    /// @brief used for convert string to dbl parameter's index
    /// @param name
    /// @return an index that denotes parameter name
    [[nodiscard]] idx_type dbl_par_str2idx(const std::string &name) const override {
        return (name == "kMethod" ? kMethod : CDASolver::dbl_par_str2idx(name));
    }

    /*--------------------------------------------------------------------------*/

    /// @brief used for convert int index idx to phrasal rapresentation
    /// @param idx
    /// @return a string that denotes index idx parameter
    [[nodiscard]] const std::string &int_par_idx2str(idx_type idx) const override {
        if (idx > intLastParLEMON_CS) {
            throw std::invalid_argument(std::to_string(idx));
        }

        static const std::string par = "kMethod";
        switch (idx) {
        case kMethod:
            return (par);
        default:
            break;
        }

        return (CDASolver::int_par_idx2str(idx));
    }

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    /// @brief used for convert dbl index idx to phrasal rapresentation
    /// @param idx
    /// @return a string that denotes index idx parameter
    [[nodiscard]] const std::string &dbl_par_idx2str(idx_type idx) const override {
        if (idx > dblLastParLEMON_CS) {
            throw std::invalid_argument(std::to_string(idx));
        }

        switch (idx) {
        default:
            break;
        }

        return (CDASolver::dbl_par_idx2str(idx));
    }

    /*--------------------------------------------------------------------------*/
    /*-------------------------------- FRIENDS
     * ---------------------------------*/
    /*--------------------------------------------------------------------------*/

    friend class MCFLemonState; // make MCFSolverState friend

    /** @}
     * ---------------------------------------------------------------------*/
    /*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
    /*--------------------------------------------------------------------------*/

  protected:
    /*--------------------------------------------------------------------------*/
    /*-------------------------- PROTECTED METHODS
     * -----------------------------*/
    /*--------------------------------------------------------------------------*/

    void guts_of_compute() override { status = f_algo->run(CSMethod(f_method)); }

    /*--------------------------------------------------------------------------*/
    /*---------------------------- PROTECTED FIELDS ---------------------------*/
    /*--------------------------------------------------------------------------*/

    std::string f_dmx_file;
    // string for DMX file output

    /*--------------------------------------------------------------------------*/
    /*--------------------- PRIVATE PART OF THE CLASS
     * --------------------------*/
    /*--------------------------------------------------------------------------*/

  private:
    /*--------------------------------------------------------------------------*/
    /*-------------------------- PRIVATE METHODS -------------------------------*/
    /*--------------------------------------------------------------------------*/

    SMSpp_insert_in_factory_h;

    /*--------------------------------------------------------------------------*/
    /*-------------------------- PRIVATE FIELDS -------------------------------*/
    /*--------------------------------------------------------------------------*/
    CSMethod f_method;
    double ticks; // Elapsed time in ticks for compute() method

    /*--------------------------------------------------------------------------*/

}; // end( class MCFLemonSolver<SMSppCostScaling, GR, V, C>  Specialization)

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS MCFLemonState ---------------------------*/
/*--------------------------------------------------------------------------*/
/// class to describe the "internal state" of a MCFLemonSolver
/** Derived class from State to describe the "internal state" of a
 * MCFLemonSolver.  */

class MCFLemonState : public State {
    /*----------------------- PUBLIC PART OF THE CLASS
     * -------------------------*/

  public:
=======
>>>>>>> 3587a0f95544d63e767d768c535e345b48386907
    /*------------- CONSTRUCTING AND DESTRUCTING MCFLemonState ----------------*/

    /// constructor, doing everything or nothing.
    /** Constructor of MCFSolverState. If provided with a pointer to a MCFSolver
     * it immediately copies its "internal state", which is the only way in
     * which the MCFSolverState can be initialised out of an existing MCFSolver.
     * If nullptr is passed (as by default), then an "empty" MCFSolverState is
     * constructed that can only be filled by calling deserialize().
     *
     * Note: to avoid having to duplicate the SMSpp_insert_in_factory_cpp call
     *       for every MCFClass, the pointer is directly that of a MCFClass,
     *       since every MCFSolver derives from a :MCFClass and we only need
     *       access to MCFGetState(). */

    /*MCFSolverState(MCFClass *mcfc = nullptr) : State()
    {
      f_state = mcfc ? mcfc->MCFGetState() : nullptr;
    }*/

    /*--------------------------------------------------------------------------*/
    /// de-serialize a MCFSolverState out of netCDF::NcGroup
    /** Should de-serialize a MCFSolverState out of netCDF::NcGroup, but in
     * fact it does not work. */

    /*void deserialize(const netCDF::NcGroup &group) override
    {
      f_state = nullptr;
    }*/

    /*--------------------------------------------------------------------------*/
    /// destructor

<<<<<<< HEAD
    // virtual ~MCFLemonState() { delete f_state; }
=======
    //virtual ~MCFLemonState() { delete f_state; }
>>>>>>> 3587a0f95544d63e767d768c535e345b48386907

    /*---------- METHODS DESCRIBING THE BEHAVIOR OF A MCFSolverState
     * -----------*/

    /// serialize a MCFSolverState into a netCDF::NcGroup
    /** The method should serialize the MCFSolverState into the provided
     * netCDF::NcGroup, so that it can later be read back by deserialize(), but
     * in fact it does not work.*/

    void serialize(netCDF::NcGroup &group) const override {}

    /*-------------------------------- FRIENDS
     * ---------------------------------*/

    template <class MCFC> friend class MCFSolver; // make MCFSolver friend

    /*-------------------- PROTECTED PART OF THE CLASS
     * -------------------------*/

  protected:
    /*-------------------------- PROTECTED METHODS
     * -----------------------------*/

    void print(std::ostream &output) const override { output << "MCFSolverState [" << this << "]"; }

    /*--------------------------- PROTECTED FIELDS
     * -----------------------------*/

<<<<<<< HEAD
    // MCFClass::MCFStatePtr f_state; ///< the (pointer to) MCFState
=======
   // MCFClass::MCFStatePtr f_state; ///< the (pointer to) MCFState
>>>>>>> 3587a0f95544d63e767d768c535e345b48386907

    /*---------------------- PRIVATE PART OF THE CLASS
     * -------------------------*/

  private:
    /*---------------------------- PRIVATE FIELDS
     * ------------------------------*/

    SMSpp_insert_in_factory_h;

    /*--------------------------------------------------------------------------*/

}; // end( class( MCFSolverState ) )

<<<<<<< HEAD
/** @} end( group( MCFSolver_CLASSES ) ) */
/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

} // namespace SMSpp_di_unipi_it
=======
  /** @} end( group( MCFSolver_CLASSES ) ) */
  /*--------------------------------------------------------------------------*/
  /*------------------- inline methods implementation ------------------------*/
  /*--------------------------------------------------------------------------*/
  //TODO: change MCFC function to Algo function.
  /*template <class Algo>
  State *MCFSolver<Algo>::get_State(void) const
  {
    return (new MCFSolverState(const_cast<MCFSolver<MCFC> *>(this)));
  }*/

  /*--------------------------------------------------------------------------*/
  //TODO: change MCFC function to Algo function.
  /*template <class Algo>
  void MCFSolver<Algo>::put_State(const State &state)
  {
    // if state is not a const MCFSolverState &, exception will be thrown
    auto s = dynamic_cast<const MCFSolverState &>(state);

    this->MCFPutState(s.f_state);
  }*/

  /*--------------------------------------------------------------------------*/
  //TODO: change MCFC function to Algo function.
  /*template <class Algo>
  void MCFSolver<Algo>::put_State(State &&state)
  {
    // if state is not a MCFSolverState &&, exception will be thrown
    auto s = dynamic_cast<MCFSolverState &&>(state);

    this->MCFPutState(s.f_state);
  }*/

  /*--------------------------------------------------------------------------*/
  //TODO: change MCFC function to Algo function.
  /*template <class Algo>
  void MCFSolver<Algo>::process_outstanding_Modification(void)
  {
    // no-frills loop: do them in order, with no attempt at optimizing
    // note that NBModification have already been dealt with and therefore need
    // not be considered here

    for (;;)
    {
      auto mod = pop();
      if (!mod)
        break;

      guts_of_poM(mod.get());
    }
  }*/ // end( MCFSolver::process_outstanding_Modification )

  /*--------------------------------------------------------------------------*/

  /*template <class Algo>
  void MCFSolver<Algo>::guts_of_poM(c_p_Mod mod)
  {
    auto MCFB = static_cast<MCFBlock *>(f_Block);
    */
    // process Modification - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /* This requires to patiently sift through the possible Modification types
     * to find what this Modification exactly is, and call the appropriate
     * method of MCFClass. */

    // GroupModification- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*if (auto tmod = dynamic_cast<const GroupModification *>(mod))
    {
      for (const auto &submod : tmod->sub_Modifications())
        guts_of_poM(submod.get());

      return;
    }*/

    // MCFBlockRngdMod- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /* Note: in the following we can assume that C, B and U are nonempty. This
     * is because they can be empty only if they are so when the object is
     * loaded. But if a Modification has been issued they are no longer empty (a
     * Modification changin nothing from the "empty" state is not issued). */

    /*if (auto tmod = dynamic_cast<const MCFBlockRngdMod *>(mod))
    {
      auto rng = tmod->rng();

      switch (tmod->type())
      {
      case (MCFBlockMod::eChgCost):
        if (rng.second == rng.first + 1)
        {
          if (!MCFB->is_deleted(rng.first))
            //TODO: change MCFC function to Algo function.
            MCFC::ChgCost(rng.first, MCFB->get_C(rng.first));
        }
        else
        {
          if (std::any_of(MCFB->get_C().data() + rng.first,
                          MCFB->get_C().data() + rng.second,
                          [](auto ci)
                          { return (std::isnan(ci)); }))
          {
            MCFBlock::Vec_CNumber NCost(MCFB->get_C().data() + rng.first,
                                        MCFB->get_C().data() + rng.second);
            for (auto &ci : NCost)
              if (std::isnan(ci))
                ci = 0;
            //TODO: change MCFC function to Algo function.
            MCFC::ChgCosts(NCost.data(), nullptr, rng.first, rng.second);
          }
          else
            //TODO: change MCFC function to Algo function.
            MCFC::ChgCosts(MCFB->get_C().data() + rng.first, nullptr,
                           rng.first, rng.second);
        }
        return;

      case (MCFBlockMod::eChgCaps):
        if (rng.second == rng.first + 1)
        {
          if (!MCFB->is_deleted(rng.first))
            //TODO: change MCFC function to Algo function.
            MCFC::ChgUCap(rng.first, MCFB->get_U(rng.first));
        }
        else
          //TODO: change MCFC function to Algo function.
          MCFC::ChgUCaps(MCFB->get_U().data() + rng.first, nullptr,
                         rng.first, rng.second);
        return;

      case (MCFBlockMod::eChgDfct):
        if (rng.second == rng.first + 1)
          //TODO: change MCFC function to Algo function.
          MCFC::ChgDfct(rng.first, MCFB->get_B(rng.first));
        else
          //TODO: change MCFC function to Algo function.
          MCFC::ChgDfcts(MCFB->get_B().data() + rng.first, nullptr,
                         rng.first, rng.second);
        return;
      //TODO: change MCFC function to Algo function.
      case (MCFBlockMod::eOpenArc):
        for (; rng.first < rng.second; ++rng.first)
          if ((!MCFB->is_deleted(rng.first)) &&
              (!MCFC::IsDeletedArc(rng.first)))
            MCFC::OpenArc(rng.first);
        return;
      //TODO: change MCFC function to Algo function.
      case (MCFBlockMod::eCloseArc):
        for (; rng.first < rng.second; ++rng.first)
          if ((!MCFB->is_deleted(rng.first)) &&
              (!MCFC::IsDeletedArc(rng.first)))
            MCFC::CloseArc(rng.first);
        return;

      case (MCFBlockMod::eAddArc):
      {
        auto ca = MCFB->get_C(rng.first);
        //TODO: change MCFC function to Algo function.
        auto arc = MCFC::AddArc(MCFB->get_SN(rng.first),
                                MCFB->get_EN(rng.first),
                                MCFB->get_U(rng.first),
                                std::isnan(ca) ? 0 : ca);
        if (arc != rng.first)
          throw(std::logic_error("name mismatch in AddArc()"));
        return;
      }

      case (MCFBlockMod::eRmvArc):
        //TODO: change MCFC function to Algo function.
        MCFC::DelArc(rng.second - 1);
        return;

      default:
        throw(std::invalid_argument("unknown MCFBlockRngdMod type"));
      }
    }

    // MCFBlockSbstMod- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (auto tmod = dynamic_cast<const MCFBlockSbstMod *>(mod))
    {
      switch (tmod->type())
      {
      case (MCFBlockMod::eOpenArc):
        for (auto arc : tmod->nms())
          if ((!MCFB->is_deleted(arc)) &&
              (!MCFC::IsDeletedArc(arc)))
            //TODO: change MCFC function to Algo function.
            MCFC::OpenArc(arc);
        return;

      case (MCFBlockMod::eCloseArc):
        for (auto arc : tmod->nms())
          if ((!MCFB->is_deleted(arc)) &&
              (!MCFC::IsDeletedArc(arc)))
            //TODO: change MCFC function to Algo function.
            MCFC::CloseArc(arc);
        return;
      }

      // have to InINF-terminate the vector of indices (damn!)
      // meanwhile, when appropriate remove the indices of deleted
      // arcs for which the operations make no sense;
      ;
      switch (tmod->type())
      {
      case (MCFBlockMod::eChgCost):
      {
        MCFBlock::Subset nmsI;
        nmsI.reserve(tmod->nms().size() + 1);
        MCFBlock::Vec_CNumber NCost;
        NCost.reserve(tmod->nms().size());
        auto &C = MCFB->get_C();
        for (auto i : tmod->nms())
          if (auto ci = C[i]; !std::isnan(ci))
          {
            NCost.push_back(ci);
            nmsI.push_back(i);
          }
        nmsI.push_back(Inf<MCFBlock::Index>());
        //TODO: change MCFC function to Algo function.
        MCFC::ChgCosts(NCost.data(), nmsI.data());
        return;
      }

      case (MCFBlockMod::eChgCaps):
      {
        MCFBlock::Subset nmsI;
        nmsI.reserve(tmod->nms().size() + 1);
        MCFBlock::Vec_FNumber NCap;
        NCap.reserve(tmod->nms().size());
        auto &C = MCFB->get_C();
        auto &U = MCFB->get_U();
        for (auto i : tmod->nms())
          if (!std::isnan(C[i]))
          {
            NCap.push_back(U[i]);
            nmsI.push_back(i);
          }
        nmsI.push_back(Inf<MCFBlock::Index>());
        //TODO: change MCFC function to Algo function.
        MCFC::ChgUCaps(NCap.data(), nmsI.data());
        return;
      }

      case (MCFBlockMod::eChgDfct):
      {
        MCFBlock::Vec_FNumber NDfct(tmod->nms().size());
        MCFBlock::Subset nmsI(tmod->nms().size() + 1);
        *copy(tmod->nms().begin(), tmod->nms().end(), nmsI.begin()) =
            Inf<MCFBlock::Index>();
        auto B = MCFB->get_B();
        for (MCFBlock::Index i = 0; i < NDfct.size(); i++)
          NDfct[i] = B[nmsI[i]];
        //TODO: change MCFC function to Algo function.
        MCFC::ChgDfcts(NDfct.data(), nmsI.data());
        return;
      }

      default:
        throw(std::invalid_argument("unknown MCFBlockSbstMod type"));
      }
    }

    // any remaining Modification is plainly ignored, since it must be an
    // "abstract" Modification, which this Solver does not need to look at

  }*/ // end( guts_of_poM )

  /*--------------------------------------------------------------------------*/
  /*--------------------------------------------------------------------------*/

} // end( namespace SMSpp_di_unipi_it )
>>>>>>> 3587a0f95544d63e767d768c535e345b48386907

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* MCFLemonSolver.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File MCFLemonSolver.h
 * ---------------------------*/
/*--------------------------------------------------------------------------*/


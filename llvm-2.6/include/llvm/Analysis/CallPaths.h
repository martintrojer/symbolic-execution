#ifndef LLVM_ANALYSIS_CALLPATHS_H
#define LLVM_ANALYSIS_CALLPATHS_H

#include <vector>
#include "llvm/Pass.h"
#include "llvm/BasicBlock.h"
#include "llvm/Function.h"
#include "llvm/Analysis/CallGraph.h"

#include <boost/config.hpp>
#include <boost/utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/depth_first_search.hpp>

using namespace boost;

namespace llvm {

  // -------(BGL)---------------------------
  //  typedef adjacency_list<setS, vecS, directedS, no_property, property<edge_weight_t, int> > Graph;
  typedef adjacency_list<setS, vecS, bidirectionalS, no_property, property<edge_weight_t, int> > Graph;
  typedef graph_traits<Graph>::vertex_descriptor Vertex;
  typedef graph_traits<Graph>::edge_descriptor Edge;
  typedef color_traits<default_color_type> Color;
  typedef std::vector<default_color_type> ColorVec;

  class my_dfs_visitor;
  
  ModulePass *createCallPathsPass(std::vector<std::vector<BasicBlock*> > *_bbpaths, std::string _filename);
  
  class CallPaths : public ModulePass {
	friend class my_dfs_visitor;
  public:
	// Do not create (createCallPathsPass is the object Factory)
	CallPaths();
	virtual void getAnalysisUsage(AnalysisUsage &AU) const;
	virtual bool runOnModule(Module &_M);

  public:
	// createCallPathsPass Input : All BB paths found
	std::vector<std::vector<BasicBlock*> > *bbpaths;
	// createCallPathsPass Input : filename to defect xml file
	std::string defectFile;
	
	// Stored pointer to the module (gotten fron runOnModule)
	Module *M;
	// LLVM Pass ID
	static char ID;

  private:
	// Set of filenames and line numbers pointing to areas of interest
	typedef std::map<std::string, std::vector<int> > defectList;
	defectList dl;
	// Parse defectFile and build the defectList (dl) map
	void getDefectList(std::string docname, defectList *res);

	//Return a Function* given a filename/srcline - used by DFS funcGraph
	Function *getFunction(std::string srcFile, int srcLine);
	//Return a BB* given a filename/srcline -- used by DFS bbGraph 
	BasicBlock *getBB(std::string srcFile, int srcLine);

	// Checks wether a BB contains the code found in a file/line compination
	bool findLineInBB(BasicBlock *BB, std::string srcFile, int srcLine);
	// Checks wether a Function contains the code found in a file/line compination
	// (uses findLineInBB)
	bool findLineInFunction(Function *F, BasicBlock **BB, std::string srcFile, int srcLine);
		
	// -------(BGL)---------------------------	
  private:
	Graph funcG, bbG;
	std::map<Function*, Vertex> funcMap;   // Map functions to vertices  
	std::map<BasicBlock*, Vertex> bbMap;

  private:
	//Return a BB* given a Vertex -- used by Dijkstra
	BasicBlock *getBB(Vertex v);
	// Given a BB, add edges to all it's successors 
	void addBBEdges(BasicBlock *BB);
	// Reset the entire weightmap.
	// FIX; no need for the BuildGraph to worry about the weightmap?
	void resetWeightMap(void);

	// Get the name of Function given a vertex (relies on the funcMap)
	std::string getName(Vertex v);
	// Get the name of BB given a vertex (relies on the bbMap)
	std::string getBBName(Vertex v);

	// Given a CallGraph (result of the CallGraph LLVM Pass) builds the
	// funcG and bbG
	void buildGraph(CallGraph *CG);
	// check wether a given paths already exists in *paths
	bool duplicate(std::vector<std::vector<Vertex> >*paths, std::vector<Vertex> p);
	// Runs Dijkstra's algortihm and find the shortest path between root and target
	void findSinglePath(std::vector<Vertex> *path, 
			    Vertex root,
			    Vertex target);
	// Tries to find a maximum of 'number' unique paths and put them in *paths
	// Runs multiple Dijsktra's and tinkers with the weightmap to find new paths
	void findPaths(std::vector<std::vector<Vertex> >*paths,
		       BasicBlock *rootBB,
		       BasicBlock *targetBB,
		       int number);

	// Helper functor used to write nice labels for functions in GraphViz dot files
  private:
	struct my_func_label_writer {
	  CallPaths *CP;
	my_func_label_writer(CallPaths *_CP) : CP(_CP) { }
	  template <class VertexOrEdge>
	  void operator()(std::ostream& out, const VertexOrEdge& v) const {
		out << "[label=\"" << CP->getName(v) << "\"]";	
	  }
	};	

	// Helper functor used to write nice labels for BBs in GraphViz dot files
	struct my_bb_label_writer {
	  CallPaths *CP;
	my_bb_label_writer(CallPaths *_CP) : CP(_CP) { }
	  template <class VertexOrEdge>
	  void operator()(std::ostream& out, const VertexOrEdge& v) const {
		out << "[label=\"" << CP->getBBName(v) << "\"]";	
	  }
	};	

  	// --- Obsolete -----
  private:
	// This is the old loop to try generating paths, doesn't work properly
	bool findBBPath(CallGraphNode *n, std::vector<BasicBlock*> &path, std::string srcFile, int srcLine);
	// BGL; find out if a BB is a pred of another (gets stuck in loops)
	bool isBBPred(BasicBlock *BB, BasicBlock *succBB);
  };
  
  // DFS Visitor class (used by the brute force DFS path builders)
  class my_dfs_visitor:public default_dfs_visitor {
  private:
	CallPaths *CP;
	Vertex target, root;
	std::vector<std::vector<Vertex> >*paths;
	std::vector<Vertex> tpath;
	std::vector<default_color_type> &colmap;
	long long ctr;
  public:
  my_dfs_visitor(CallPaths *_CP, Vertex _target, Vertex _root,
		 std::vector<std::vector<Vertex> >*_paths,
		 std::vector<default_color_type> &_colmap) :
	CP(_CP), target(_target), root(_root), paths(_paths), colmap(_colmap), ctr(0) { }
	bool duplicate(std::vector<Vertex> p)
	{
	  for (std::vector<std::vector<Vertex> >::iterator it=paths->begin(); it != paths->end(); ++it) {
		std::vector<Vertex> cpath = *it;
		if (cpath.size() != p.size())
		  continue;
		unsigned int i;
		for (i=0; i<cpath.size(); ++i) {
		  if (p[i] != cpath[i])
			break;
		}
		if (i==cpath.size())
		  return true;
	  }
	  return false;
	}
	void newPath(void)
	{
	  if (tpath.size() > 0)
		if (tpath.front() == root)
		  if (!duplicate(tpath))
			paths->push_back(tpath);

	  std::cerr << "newPath " << paths->size() << " : " << tpath.size() << "\n";
	}
	void discover_vertex(Vertex u, const Graph & g)
	{
	  //	  std::cerr << "discover " << CP->getName(u) << "\n";

	  ctr++;

	  if ((ctr%1000000)==0)
	    std::cerr << "visited " << ctr << " vertices -- tpath size " 
		      << tpath.size() << "\n";

	  tpath.push_back(u);
	  if (u == target) {
		newPath();
	  }
	}
	void finish_vertex(Vertex u, const Graph & g)
	{
	  //	  std::cerr << "finish " << CP->getName(u) << " : " << colmap[u] <<"\n";

	  tpath.pop_back();
	  colmap[u] = Color::white();	  
	}
  };
}

#endif

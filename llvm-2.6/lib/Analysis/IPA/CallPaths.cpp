// export CXXFLAGS="-I /usr/include/libxml2 -fexceptions"
//
// Number of Paths found when run from klee is somewhat bizarre.
// Combinations of the "--libc=uclibc --posix-runtime" and "--optimize" effects it,
// but not really as you think; with optimize often gives more paths.

#include <iostream>
#include <fstream>

#include "llvm/Analysis/CallPaths.h"

#include "llvm/Module.h"
#include "llvm/Support/GraphWriter.h"
#include "llvm/Pass.h"
#include "llvm/Value.h"
#include "llvm/Analysis/CallGraph.h"
#include "llvm/Support/CFG.h"
#include "llvm/Analysis/Dominators.h"
#include "llvm/IntrinsicInst.h"
#include "llvm/Analysis/ValueTracking.h"

using namespace llvm;

// This will add -lxml2 dependancy when linking, and will break LLVM builds
// use CXXFLAGS="-I /usr/include/libxml2" and LDFLAGS "-l xml2"
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>

// Uses Boost Graph Library
// FIX: There is something really ugly going on here...
// In order to being able to link with libLLVMipa.a we have to hack
// /usr/include/boost/serialization/throw_exception.hpp
#include <boost/config.hpp>
#include <boost/utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

using namespace boost;

typedef std::map<std::string, std::vector<int> > defectList;

// -----------------------------------------------------------
// --- LLVM Pass Stuff

ModulePass *llvm::createCallPathsPass(std::vector<std::vector<BasicBlock*> > *_bbpaths, std::string _filename)
{
  // One of the few news.. Trying to keep everything on the stack
  CallPaths *cg = new CallPaths();
  cg->bbpaths = _bbpaths;
  cg->defectFile = _filename;
  return cg;
}

char CallPaths::ID = 0;
static RegisterPass<CallPaths>
X("callgraph-cfg", "Do CallGraph stuff", false, true);

CallPaths::CallPaths() : ModulePass(&ID), bbpaths(NULL) { }

void CallPaths::getAnalysisUsage(AnalysisUsage &AU) const {
  AU.addRequired<CallGraph>();
  AU.setPreservesAll();
}

bool CallPaths::runOnModule(Module &_M) {
  M = &_M;
  CallGraph *CG = &getAnalysis<CallGraph>();
  CallGraphNode *root = CG->getRoot();
  if (root == NULL) return false;
  Function *rootF = root->getFunction();
  if (rootF == NULL) return false;
  BasicBlock *rootBB = &rootF->getEntryBlock();
  if (rootBB == NULL) return false;

  // Parse and extract data from Defect XML
  getDefectList(defectFile, &dl);
  if (dl.size()==0) return true;

  // Build the BGL Graph
  buildGraph(CG);

  // --------------------------------------------------------------
  // ---- Brute Force DFS on bbGraph 
#if 0
  Vertex s = bbMap[rootBB];
  for (defectList::iterator dit=dl.begin(); dit != dl.end(); ++dit) {
	std::string file = dit->first;
	std::vector<int> lines = dit->second;
	for (std::vector<int>::iterator lit=lines.begin(); lit != lines.end(); ++lit) {
	  std::cerr << "Looking for '" << file << "' (" << *lit << ")\n";
	  BasicBlock *BB = getBB(file, *lit);
	  if (BB!= NULL) {
		Function *F = BB->getParent();
	    std::cerr << "Finding Paths to " << F->getNameStr() << " : " << BB->getNameStr() << "\n";
	    std::vector<default_color_type> colmap(num_vertices(bbG));
	    std::vector<std::vector<Vertex> >paths;
	    my_dfs_visitor vis(this, bbMap[BB], s, &paths, colmap);
	    depth_first_search(bbG, vis, &colmap[0], s);

		// Now convert to bbPaths...
		
	    for (unsigned int i = 0; i < paths.size(); ++i) {
		  std::vector<BasicBlock*> tpath;
	      for (std::vector<Vertex>::iterator it=paths[i].begin(); it != paths[i].end(); ++it) {
			//			std::cout << getBBName(*it) << " -- ";
			tpath.push_back(getBB(*it));
	      }
		  bbpaths->push_back(tpath);
		  //	      std::cout << "\n";
	    }
	  }
	}
  }
#endif

  // --------------------------------------------------------------
  // ---- Brute Force DFS on funcGraph 
#if 0
  Vertex s = funcMap[rootF];
  for (defectList::iterator dit=dl.begin(); dit != dl.end(); ++dit) {
	std::string file = dit->first;
	std::vector<int> lines = dit->second;
	for (std::vector<int>::iterator lit=lines.begin(); lit != lines.end(); ++lit) {
	  std::cerr << "Looking for '" << file << "' (" << *lit << ")\n";
	  Function *F = getFunction(file, *lit);
	  if (F!= NULL) {
	    std::cerr << "Finding Paths to " << F->getNameStr() << "\n";
	    std::vector<default_color_type> colmap(num_vertices(funcG));
	    std::vector<std::vector<Vertex> >paths;
	    my_dfs_visitor vis(this, funcMap[F], s, &paths, colmap);
	    depth_first_search(funcG, vis, &colmap[0], s);

	    for (unsigned int i = 0; i < paths.size(); ++i) {
	      for (std::vector<Vertex>::iterator it=paths[i].begin(); it != paths[i].end(); ++it) {
			std::cout << getName(*it) << " -- ";
	      }
	      std::cout << "\n";
	    }
	  }
	}
  }
#endif

  // --------------------------------------------------------------
  // ---- "Dijkstra solution"
  // Requires "-fexceptions" in CXXFLAGS
#if 1

  for (defectList::iterator dit=dl.begin(); dit != dl.end(); ++dit) {
    std::string file = dit->first;
    std::vector<int> lines = dit->second;
    for (std::vector<int>::iterator lit=lines.begin(); lit != lines.end(); ++lit) {
      std::cerr << "Looking for '" << file << "' (" << *lit << ")\n";
      BasicBlock *BB = getBB(file, *lit);
      if (BB==NULL)
		continue;

      // Now find a shitload of paths for this targetBB

      std::vector<std::vector<Vertex> > paths;
      findPaths(&paths, rootBB, BB, 1000);

      // Convert to BB lists
      for (std::vector<std::vector<Vertex> >::iterator pit=paths.begin();
		   pit != paths.end(); ++pit) {
		std::vector<BasicBlock *> tpath;
		for (std::vector<Vertex>::iterator vit=pit->begin();
			 vit != pit->end(); ++vit) {
		  tpath.push_back(getBB(*vit));
		}	
		bbpaths->push_back(tpath);
      }
    }
  }

  // FIX: We might want to add a second pass here adding "loops" into the paths?
  // maybe strongly_connected algorithm can be used for this?
  // Well, if you have loops you no longer have a path, but another graph,
  // which makes the Searcher more complicated... but maybe required?
#endif

  std::cerr << "Found " << bbpaths->size() << " path(s)\n";
  return false;
}

// -----------------------------------------------------------
// --- XML parsing funciton (relies on libxml2)

void CallPaths::getDefectList(std::string docname, defectList *res)
{
  xmlDocPtr doc;
  xmlNodePtr cur;

  doc = xmlParseFile(docname.c_str());
	
  if (!doc) return;
	
  cur = xmlDocGetRootElement(doc);
	
  if (!cur) {
	xmlFreeDoc(doc);
	return;
  }

  // Iterate over defects
  while (cur) {
	if (!xmlStrcmp(cur->name, (const xmlChar *) "defect")) {
	  xmlNodePtr d = cur->xmlChildrenNode;
	  std::string file;
	  std::vector<int> lines;
	  while (d) {
		if (!xmlStrcmp(d->name, (const xmlChar *) "file")) {
		  file = (char*)xmlNodeListGetString(doc, d->xmlChildrenNode, 1);
		  lines.clear();
		}
		else if (!xmlStrcmp(d->name, (const xmlChar *) "event")) {
		  xmlNodePtr e = d->xmlChildrenNode;
		  while (e) {
			if (!xmlStrcmp(e->name, (const xmlChar *) "line")) {
			  char *val = (char*)xmlNodeListGetString(doc, e->xmlChildrenNode, 1);
			  lines.push_back(atoi(val));
			}
			e = e->next;
		  }
		}
		d = d->next;
	  }
	  if (file != "" && lines.size() > 0)
		res->insert(std::make_pair(file,lines));
	}
	cur = cur->next;
  }
	
  xmlFreeDoc(doc);
}

// -----------------------------------------------------------
// --- LLVM Function/BB help functions

// Converts a DSPI to a filename with PATH (requires debug data in the .o file)
static std::string getDSPIPath(DbgStopPointInst *dspi)
{
  std::string dir, file;
  bool res = GetConstantStringInfo(dspi->getDirectory(), dir);
  assert(res && "GetConstantStringInfo failed");
  res = GetConstantStringInfo(dspi->getFileName(), file);
  assert(res && "GetConstantStringInfo failed");
  if (dir.empty()) {
	return file;
  } else if (*dir.rbegin() == '/') {
	return dir + file;
  } else {
	return dir + "/" + file;
  }
}

bool CallPaths::findLineInBB(BasicBlock *BB, std::string srcFile, int srcLine)
{
  for (BasicBlock::iterator it = BB->begin(), ie = BB->end(); it != ie; ++it) {
	if (DbgStopPointInst *dspi = dyn_cast<DbgStopPointInst>(&*it)) {
	  std::string bbFile = getDSPIPath(dspi);
	  int bbLine = dspi->getLine();
	  //	  std::cerr << " :checking " << bbFile << " : " << bbLine << "\n";
	  if ((bbFile == srcFile) && (bbLine == srcLine))
		return true;
	}
  }		
  return false;
}

bool CallPaths::findLineInFunction(Function *F, BasicBlock **BB, std::string srcFile, int srcLine)
{
  for (Function::iterator bbIt = F->begin(), bb_ie = F->end(); bbIt != bb_ie; ++bbIt) {
	*BB = bbIt;
	if (findLineInBB(*BB,srcFile,srcLine))
	  return true;
  }
  return false;
}

Function *CallPaths::getFunction(std::string srcFile, int srcLine)
{
  for (Module::iterator fit=M->begin(); fit != M->end(); ++fit) {
    BasicBlock *BB;
    if (findLineInFunction(fit, &BB, srcFile, srcLine))
      return fit;
  }
  return NULL;
}

BasicBlock *CallPaths::getBB(std::string srcFile, int srcLine)
{
  for (Module::iterator fit=M->begin(); fit != M->end(); ++fit) {
	Function *F = fit;
	for (Function::iterator bbit=F->begin(); bbit != F->end(); ++bbit) {
	  BasicBlock *BB = bbit;
	  if (findLineInBB(BB, srcFile, srcLine))
		return BB;
	}
  }
  return NULL;
}

// -----------------------------------------------------------
// --- BGL Help functions

BasicBlock *CallPaths::getBB(Vertex v)
{
  for (std::map<BasicBlock*, Vertex>::iterator it=bbMap.begin(); it != bbMap.end(); ++it) {
	if (v == it->second)
	  return it->first;
  } 
  return NULL;
}

// This functions relies on the fact that a OutEdgeList of type setS doesn't allow
// non-unique edges (and returns false if such an edge is being inserted)
void CallPaths::addBBEdges(BasicBlock *BB)
{
  graph_traits<Graph>::edge_descriptor e; bool inserted;
  property_map<Graph, edge_weight_t>::type bbWeightmap = get(edge_weight, bbG);

  for (succ_iterator si = succ_begin(BB); si != succ_end(BB); ++si) {
	//	const char *bbName = si->getNameStr().c_str();
	// Dont get stuck in loops
    // if (!isBBPred(BB,*si))
	boost::tie(e, inserted) = add_edge(bbMap[BB], bbMap[*si], bbG);
	if (inserted)
	  addBBEdges(*si);
	bbWeightmap[e] = 1;
  }
}

void CallPaths::resetWeightMap(void)
{
  property_map<Graph, edge_weight_t>::type funcWeightmap = get(edge_weight, funcG);
  property_map<Graph, edge_weight_t>::type bbWeightmap = get(edge_weight, bbG);
  graph_traits<Graph>::edge_iterator ei, ei_end;

#if 0
  for (boost::tie(ei,ei_end) = edges(funcG); ei != ei_end; ++ei)
    funcWeightmap[*ei] = 1;
#endif

  for (boost::tie(ei,ei_end) = edges(bbG); ei != ei_end; ++ei)
    bbWeightmap[*ei] = 1;
}

std::string CallPaths::getName(Vertex v)
{
  std::string res = "<null>";
  
  for (std::map<Function *, Vertex>::iterator it=funcMap.begin(); it != funcMap.end(); ++it) {
	if (v == it->second) {
	  Function *F = it->first;
	  if (F != NULL)
		res = it->first->getNameStr();
	}
  } 
  return res;
}

std::string CallPaths::getBBName(Vertex v)
{
  std::string res = "<null>";
  
  for (std::map<BasicBlock*, Vertex>::iterator it=bbMap.begin(); it != bbMap.end(); ++it) {
	if (v == it->second) {
	  BasicBlock *BB = it->first;
	  if (BB != NULL) {
	    Function *F = BB->getParent();
	    res = F->getNameStr() + ":" + BB->getNameStr();
	  }
	}
  } 
  return res;
}

void CallPaths::buildGraph(CallGraph *CG)
{
  std::cerr << "Building Vertices... ";
  for (Module::iterator fit=M->begin(); fit != M->end(); ++fit) {
    Function *F = fit;
    funcMap[F] = add_vertex(funcG);
    for (Function::iterator bbIt = F->begin(), bb_ie = F->end(); bbIt != bb_ie; ++bbIt) {
      BasicBlock *BB = bbIt;
      bbMap[BB] = add_vertex(bbG);
    }
  }    
  std::cerr << "funcMap: " << num_vertices(funcG) << " - bbMap: " << num_vertices(bbG) << "\n";

  property_map<Graph, edge_weight_t>::type funcWeightmap = get(edge_weight, funcG);
  property_map<Graph, edge_weight_t>::type bbWeightmap = get(edge_weight, bbG);
  std::cerr << "Building Edges... ";

  for (Module::iterator fit = M->begin(); fit != M->end(); ++fit) {
    Function *F = fit;
	//	const char *fName = F->getNameStr().c_str();
    CallGraphNode *cgn = CG->getOrInsertFunction(F);
    if (cgn == NULL)
      continue;

	graph_traits<Graph>::edge_descriptor e; bool inserted;

	// Create edges for Function-internal BBs
    if (!F->empty()) {
      BasicBlock *BB = &F->getEntryBlock();
	  addBBEdges(BB);
    }

	// Create cross-functional edges
    for (CallGraphNode::iterator cit = cgn->begin(); cit != cgn->end(); ++cit) {
      CallGraphNode *tcgn = cit->second;
      Function *tF = tcgn->getFunction();
      if (tF == NULL)
		continue;
      //	  const char *tfName = tF->getNameStr().c_str();
      boost::tie(e, inserted) = add_edge(funcMap[F], funcMap[tF], funcG);
      funcWeightmap[e] = 1;

	  if (tF->empty())
		continue;

	  CallSite myCs = cit->first;
      Instruction *myI = myCs.getInstruction();
      BasicBlock *myBB = myI->getParent();

      Function::iterator cBBit = tF->getEntryBlock();
      BasicBlock *childBB = cBBit;
      if (childBB == NULL) 
		continue;
      boost::tie(e, inserted) = add_edge(bbMap[myBB], bbMap[childBB], bbG);
      bbWeightmap[e] = 1;
    }
  }  

  std::cerr << "funcMap: " << num_edges(funcG) << " - bbMap: " << num_edges(bbG) << "\n";

  std::ofstream funcs_file("funcs.dot");
  boost::write_graphviz(funcs_file, funcG, my_func_label_writer(this));    

  //  std::ofstream bbs_file("bbs.dot");
  //  boost::write_graphviz(bbs_file, bbG, my_bb_label_writer(this));    
}

// FIX: Are we spending alot of time in these loops?
// Area for optimization...
bool CallPaths::duplicate(std::vector<std::vector<Vertex> >*paths, std::vector<Vertex> p)
{
  if (p.size() == 0)
	return true;
  for (std::vector<std::vector<Vertex> >::iterator it=paths->begin(); 
       it != paths->end(); ++it) {
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

void CallPaths::findSinglePath(std::vector<Vertex> *path, 
							   Vertex root,
							   Vertex target)
{
  std::vector<Vertex> p(num_vertices(bbG));
  std::vector<int> d(num_vertices(bbG));  
  property_map<Graph, vertex_index_t>::type indexmap = get(vertex_index, bbG);
  property_map<Graph, edge_weight_t>::type bbWeightmap = get(edge_weight, bbG);

  dijkstra_shortest_paths(bbG, root, &p[0], &d[0], bbWeightmap, indexmap, 
						  std::less<int>(), closed_plus<int>(), 
						  (std::numeric_limits<int>::max)(), 0,
						  default_dijkstra_visitor());
  
  //  std::cout << "shortest path:" << std::endl;
  while (p[target] != target) {
    //    std::cout << getBBName(v) << " - ";
    path->insert(path->begin(), target);
    target = p[target];
  }
  // Put the root in the list aswell since the loop above misses that one
  if (!path->empty())
    path->insert(path->begin(), root);
}

void CallPaths::findPaths(std::vector<std::vector<Vertex> >*paths,
						  BasicBlock *rootBB,
						  BasicBlock *targetBB,
						  int number)
{
  Vertex root = bbMap[rootBB];
  Vertex target = bbMap[targetBB];
  property_map<Graph, edge_weight_t>::type bbWeightmap = get(edge_weight, bbG);
  int ctr=0;
  int last_path=0, last_edge = 0;
  int incr = num_edges(bbG) / 4;
 
  resetWeightMap();
  bool first = true;

  while (ctr < number) {
    std::vector<Vertex> tpath;
    findSinglePath(&tpath, root, target);
    if (duplicate(paths, tpath)) {      
      
      bool done = false, next = false;

	  int pctr=0;
      // Find an edge in previous path(s) and tinker with weight
      for (std::vector<std::vector<Vertex> >::iterator pit=paths->begin();
		   (!done) && (pit != paths->end()); ++pit, pctr++) {
		int ectr=0;
		for (std::vector<Vertex>::reverse_iterator vit=pit->rbegin();
			 (!done) && (vit != pit->rend()); ++vit, ectr++) {

		  //		  std::cout << "Path: " << pctr << " Edge: " << ectr << "\n";
		  
		  // FIX: There is an iterator for this, check out inv_adjacency_iterator

		  std::vector<Vertex>::reverse_iterator predit = vit+1;
		  if (predit == pit->rend())
			continue; //No more edges in this path
		  graph_traits<Graph>::vertex_descriptor v, predv;
		  v = *vit;
		  predv = *predit;
		  bool flag;
		  graph_traits<Graph>::edge_descriptor e;
		  boost::tie(e,flag) = edge(predv, v, bbG);


		  
		  // We need to only color one edge at the time here and not
		  // incrementally all of them. This means that we need to remember next
		  // edge to color...
		  if (flag && first) {
			// first run?
			bbWeightmap[e] += incr;
			first = false;
			done = true;
		  }
		  else if (flag && (pctr==last_path) && (ectr==last_edge)) {
			// Uncomment if you want behaviour as explained above...
			//			bbWeightmap[e] = 1;
			next = true;
		  }
		  else if (flag && next) {
			//			std::cerr << "findPaths: Increasing weight of (" << predv << "," << v << ") by " << incr << "\n";
			bbWeightmap[e] += incr;
			last_path = pctr;
			last_edge = ectr;
			done = true;
		  }
		}
      }      
      if (!done) {
		std::cerr << "findPaths: Bailing out, no more edges to tinker with\n";
		return; // we are totally done now, no more edges to tinker with
      }
    }
    else {
      std::cerr << "findPaths: Found path with length " << tpath.size()
				<< " {processing " << last_path << " of " << paths->size()+1 << "}\n";
      paths->push_back(tpath);
	  incr = tpath.size()/1;    // This should be a more realistict number then the one guessed above
	                            // FIX: Tinkering with the divisor generates different sets (and number of) paths
      ctr++;
    }
  }
}

// ----------------------------------------------------------------
// --- Obsolete stuff 
// ----------------------------------------------------------------

// FIX: Remove, is not used...
// BGL; find out if a BB is a pred of another (gets stuck in loops)
bool CallPaths::isBBPred(BasicBlock *BB, BasicBlock *succBB)
{
  for (pred_iterator pi = pred_begin(BB); pi != pred_end(BB); ++pi) {
	if (succBB == *pi)
	  return true;
  }
  return false;
}

// FIX: Obsolete (and doesn't really work)
// Gets stuck in loops etc
bool CallPaths::findBBPath(CallGraphNode *n, std::vector<BasicBlock*> &path, std::string srcFile, int srcLine)
{
  if (n == NULL) return false;

  Function *F = n->getFunction();

  std::cerr << "Processing " << F->getNameStr() << "\n";

  // Are we on a leaf?
  if (n->size() == 0) {
	BasicBlock *bb=NULL;
	if (findLineInFunction(F,&bb,srcFile,srcLine)) {
	  path.push_back(bb);
	  return true;
	}
  }
  
  for (CallGraphNode::iterator it = n->begin(); it != n->end(); ++it) {
	CallSite cs = it->first;
	CallGraphNode *tCGN = it->second;
	Instruction *tI = cs.getInstruction();
	if (tI == NULL) return false;
	BasicBlock *bb = tI->getParent();
	Function *tF = tCGN->getFunction();

	path.push_back(bb);
	if (findLineInBB(bb,srcFile,srcLine))
	  return true;

	if (tF != F) {    // Dont get stuck in recursion
	  if (findBBPath(tCGN,path,srcFile,srcLine))
		return true;
	}

	std::cerr << " Dead end, reverting...\n";  // FIX: This is misleading, not really correct.
	path.pop_back();
  }
  return false;
}

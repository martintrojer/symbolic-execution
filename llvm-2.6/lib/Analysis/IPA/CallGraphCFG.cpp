
#include "llvm/Analysis/CallGraphCFG.h"

#include "llvm/Module.h"
#include "llvm/Support/GraphWriter.h"
#include "llvm/Pass.h"
#include "llvm/Value.h"
#include "llvm/Analysis/CallGraph.h"
#include "llvm/Analysis/Dominators.h"
#include <iostream>
#include <fstream>

#include "llvm/IntrinsicInst.h"
#include "llvm/Analysis/ValueTracking.h"

using namespace llvm;

typedef std::map<std::string, std::vector<int> > defectList;

ModulePass *llvm::createCallGraphCFGPass(std::vector<std::vector<BasicBlock*> > *_bbpaths, std::string _filename)
{
  CallGraphCFG *cg = new CallGraphCFG();
  cg->bbpaths = _bbpaths;
  cg->defectFile = _filename;
  return cg;
}

char CallGraphCFG::ID = 0;
static RegisterPass<CallGraphCFG>
X("callgraph-cfg", "Do CallGraph stuff", false, true);

CallGraphCFG::CallGraphCFG() : ModulePass(&ID), bbpaths(NULL) { }

void CallGraphCFG::getAnalysisUsage(AnalysisUsage &AU) const {
  AU.addRequired<CallGraph>();
  AU.setPreservesAll();
}

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

bool CallGraphCFG::findLineInBB(BasicBlock *BB, std::string srcFile, int srcLine)
{
  for (BasicBlock::iterator it = BB->begin(), ie = BB->end(); it != ie; ++it) {
	if (DbgStopPointInst *dspi = dyn_cast<DbgStopPointInst>(&*it)) {
	  std::string bbFile = getDSPIPath(dspi);
	  int bbLine = dspi->getLine();
	  std::cerr << " :checking " << bbFile << " : " << bbLine << "\n";
	  if ((bbFile == srcFile) && (bbLine == srcLine))
		return true;
	}
  }		
  return false;
}

bool CallGraphCFG::findLineInFunction(Function *F, BasicBlock **BB, std::string srcFile, int srcLine)
{
  for (Function::iterator bbIt = F->begin(), bb_ie = F->end(); bbIt != bb_ie; ++bbIt) {
	*BB = bbIt;
	if (findLineInBB(*BB,srcFile,srcLine))
		return true;
  }
  return false;
}


bool CallGraphCFG::findBBPath(CallGraphNode *n, std::vector<BasicBlock*> &path, std::string srcFile, int srcLine)
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

bool CallGraphCFG::runOnModule(Module &M) {
  CallGraph *CG = &getAnalysis<CallGraph>();
  CallGraphNode *root = CG->getRoot();
  if (root == NULL) return true;

  getDefectList(defectFile, &dl);
  if (dl.size()==0) return true;

  for (defectList::iterator dit=dl.begin(); dit != dl.end(); ++dit) {
	std::string file = dit->first;
	std::vector<int> lines = dit->second;
	for (std::vector<int>::iterator lit=lines.begin(); lit != lines.end(); ++lit) {
	  std::cerr << "Mapping defect to BBs " << file << ":" << *lit << "\n";
	  std::vector<BasicBlock*> tpath;
	  findBBPath(root, tpath, file, *lit);
	  if (tpath.size() == 0)
		std::cerr << "! unreachable " << file << ":" << *lit << "\n";
	  else
		bbpaths->push_back(tpath);
	}
  }
  return false;
}

// This will add -lxml2 dependancy when linking, and will break LLVM builds
// use CXXFLAGS "-I /usr/include/libxml2" and LDFLAGS "-l xml2"

#include <libxml/xmlmemory.h>
#include <libxml/parser.h>

void CallGraphCFG::getDefectList(std::string docname, defectList *res)
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

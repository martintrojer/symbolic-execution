
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

namespace {
  struct CallGrapCFG : public ModulePass {
    static char ID; // Pass ID, replacement for typeid
    CallGrapCFG() : ModulePass(&ID) {}

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
	
	static void printBB(Function *F)
	{
	  for (Function::iterator bbIt = F->begin(), bb_ie = F->end(); bbIt != bb_ie; ++bbIt) {
		for (BasicBlock::iterator it = bbIt->begin(), ie = bbIt->end(); it != ie; ++it) {
		  if (DbgStopPointInst *dspi = dyn_cast<DbgStopPointInst>(&*it)) {
			std::string file = getDSPIPath(dspi);
			int line = dspi->getLine();
			std::cout << "  " << bbIt->getNameStr() << " '" << file << "' (" << line << ")\n";
			break;
		  }		
		}
	  }
	}

	static void printCG(CallGraphNode *n, int depth)
	{
	  if (n == NULL) return;

	  Function *F = n->getFunction();
	  if (F == NULL) return;

	  std::string dstr;

	  for (int i=0; i<depth; i++) dstr += "-";
	  
	  std::cout << "[" << dstr << "] " << F->getNameStr() << "\n";
	  printBB(F);

	  for (CallGraphNode::const_iterator it=n->begin(); it != n->end(); ++it) {
		CallGraphNode *tCGN = it->second;
		Function *tF = tCGN->getFunction();
		if (tF != F)
		  printCG(it->second, depth+1);
	  }
	}
	
    virtual bool runOnModule(Module &M) {
	  CallGraph *CG = &getAnalysis<CallGraph>();
	  //const Function *mainFn = M.getFunction("main");

	  CallGraphNode *root = CG->getRoot();
	  if (root == NULL) return true;

	  Function *F = root->getFunction();
	  if (F == NULL) return true;

	  std::cout << "CallGraph for " << F->getNameStr() << "\n";

	  printCG(root, 0);
	  
	  if(0) // All Function Names
		{
		  for (CallGraph::const_iterator it=CG->begin(); it != CG->end(); ++it) {
			CallGraphNode *CGN = it->second;
			Function *F = CGN->getFunction();
			if (F)
			  std::cout << F->getNameStr() << "\n";
		  }
		}
	  
      return false;
    }

    void print(std::ostream &OS) const {}
    void print(std::ostream &OS, const llvm::Module*) const {}

    virtual void getAnalysisUsage(AnalysisUsage &AU) const {
      AU.addRequired<CallGraph>();
      AU.setPreservesAll();
    }
  };

  char CallGrapCFG::ID = 0;
  RegisterPass<CallGrapCFG> P2("callgraph-cfg",
							   "Do CallGraph stuff");
}


#ifndef LLVM_ANALYSIS_CALLGRAPHCFG_H
#define LLVM_ANALYSIS_CALLGRAPHCFG_H

#include <vector>
#include "llvm/Pass.h"
#include "llvm/BasicBlock.h"
#include "llvm/Function.h"
#include "llvm/Analysis/CallGraph.h"

namespace llvm {

  class CallGraphCFG : public ModulePass {
  public:
	typedef std::map<std::string, std::vector<int> > defectList;

  public:
	std::vector<std::vector<BasicBlock*> > *bbpaths;
	std::string defectFile;
	static char ID;
	defectList dl;

	CallGraphCFG();
	virtual void getAnalysisUsage(AnalysisUsage &AU) const;
	virtual bool runOnModule(Module &M);

  private:
	// At function level now, should rellay be at BB level
	bool findLineInBB(BasicBlock *BB, std::string srcFile, int srcLine);
	bool findLineInFunction(Function *F, BasicBlock **BB, std::string srcFile, int srcLine);
	bool findBBPath(CallGraphNode *n, std::vector<BasicBlock*> &path, std::string srcFile, int srcLine);

	void getDefectList(std::string docname, defectList *res);
  };

  ModulePass *createCallGraphCFGPass(std::vector<std::vector<BasicBlock*> > *_bbpaths, std::string _filename);
}

#endif

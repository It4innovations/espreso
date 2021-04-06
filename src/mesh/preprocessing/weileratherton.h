
#ifndef SRC_MESH_PREPROCESSING_WEILERATHERTON_H_
#define SRC_MESH_PREPROCESSING_WEILERATHERTON_H_

#include "basis/containers/point.h"
#include <vector>
#include <iostream>
#include <iomanip>

namespace espreso {

typedef enum {
	Left,
	On,
	Right
} LeftOnRight;

inline const char* getLeftOnRightName(LeftOnRight leftOnRight) {
	switch (leftOnRight) {
	case Left:  return "L";
	case On:    return "O";
	case Right: return "R";
	default:    return " ";
	}
}

typedef enum {
	PList,
	QList,
	CList
} NodeTypes;

inline const char* getNodeTypeName(NodeTypes nodeType) {
	switch (nodeType) {
	case PList: return "PList";
	case QList: return "QList";
	case CList: return "CList";
	default:    return "     ";
	}
}

typedef enum {
	Single,
	PDouble,
	QDouble,
	PQDouble,
	None
} CrossTypes;

inline const char* getCrossTypeName(CrossTypes crossType) {
	switch (crossType) {
	case   Single: return "  Single";
	case  PDouble: return " PDouble";
	case  QDouble: return " QDouble";
	case PQDouble: return "PQDouble";
	default:       return "        ";
	}
}

typedef enum {
	PLine,
	QLine,
	PQLine,
	PmQiLine,
	NLine
} NextLineTypes;

inline const char* getNextLineTypeName(NextLineTypes nextLine) {
	switch (nextLine) {
	case    PLine: return "   PLine";
	case    QLine: return "   QLine";
	case   PQLine: return "  PQLine";
	case PmQiLine: return "PmQiLine";
	default:       return "        ";
	}
}



struct polyItem {
	int                       wn;           // wn={0,2,..} ... is outside;   wn=1 ... is in other poly  see http://geomalgorithms.com/a03-_inclusion.html
	bool                      inside;       // from wn
	bool                      insidesetted; // from cross points
	std::vector<unsigned int> inCross;      // indices to incidence crosspoints
	friend std::ostream& operator<<(std::ostream& os, const polyItem& );
};

struct crossItem {
	double        x;
	double        y;
	CrossTypes    type;

	unsigned int  i;
	unsigned int  j;
	unsigned int  k;
	double        pPar;
	unsigned int  pNext;
	NodeTypes     pNextType;

	unsigned int  m;
	unsigned int  n;
	unsigned int  o;
	double        qPar;
	unsigned int  qNext;
	NodeTypes     qNextType;

	NextLineTypes nextLine;

	friend std::ostream& operator<<(std::ostream& os, const crossItem& );
};

struct nodeState {
	NodeTypes nodeType;       // in what node type we are
	unsigned int listInd;     // actual index to plist, qlist or clist
	friend std::ostream& operator<<(std::ostream& os, const nodeState& );
};

class WeilerAthertonState {
	std::vector< polyItem> & pList;
	std::vector< polyItem> & qList;
	std::vector<crossItem> & cList;
	std::vector<bool> pVisited;
	std::vector<bool> qVisited;
	std::vector<bool> cVisited;
	unsigned int pVisitedCountdown;
	unsigned int cVisitedCountdown;
	bool notEnd;                // if still continue
	bool processingOutput;      // if we are outputting poly
	nodeState processedNode;
	nodeState outputStartNode;
	NextLineTypes previousLineType;
public:
	WeilerAthertonState(std::vector<polyItem> &, std::vector<polyItem> &, std::vector<crossItem> &);
	inline bool getNotEnd() {if (pVisitedCountdown<=0 && cVisitedCountdown<=0 && !processingOutput) {notEnd = false;} return notEnd;}
	inline nodeState getProcessedNode()   {return   processedNode;}
	inline nodeState getOutputStartNode() {return outputStartNode;}
	inline bool isProcessingOutput() {return processingOutput;}
	inline NodeTypes getProcessedNodeType() {return processedNode.nodeType;}
	inline unsigned int getProcessedNodeInd() {return processedNode.listInd;}
	inline unsigned int getPVisitedCountdown() {return pVisitedCountdown;}
	inline unsigned int getCVisitedCountdown() {return cVisitedCountdown;}
	inline std::vector<bool> getPVisited() {return pVisited;}
	inline std::vector<bool> getCVisited() {return cVisited;}
	inline void setEnd() {notEnd = false; }
	inline bool wasVisited() {switch(processedNode.nodeType){case PList: return pVisited[processedNode.listInd];        break; case QList: return qVisited[processedNode.listInd];        break; default: return cVisited[processedNode.listInd]; break;}}
	inline bool isInside() {  switch(processedNode.nodeType){case PList: return pList[   processedNode.listInd].inside; break; case QList: return qList[   processedNode.listInd].inside; break; default: return true;                            break;}}
	void pNext();
	void qNext();
	void cNext();
	void startProcessOutput();
	bool canStartProcessOutput();
	void processOutput();
	void endProcessOutput();
	inline bool processedNodeIsOutputStartNode() { return (processedNode.nodeType == outputStartNode.nodeType) && (processedNode.listInd == outputStartNode.listInd);}
	friend std::ostream& operator<<(std::ostream& os, const WeilerAthertonState& );
};

CrossTypes    decideCrossType(double, double, double);
NextLineTypes decideSingleCrossNextLineType(  LeftOnRight, LeftOnRight);
NextLineTypes decidePDoubleCrossNextLineType( LeftOnRight, LeftOnRight, LeftOnRight, LeftOnRight, LeftOnRight, bool);
NextLineTypes decideQDoubleCrossNextLineType( LeftOnRight, LeftOnRight, LeftOnRight, bool);
NextLineTypes decidePQDoubleCrossNextLineType(LeftOnRight, LeftOnRight, LeftOnRight, LeftOnRight, LeftOnRight, LeftOnRight, LeftOnRight, LeftOnRight, bool);

void pushBackSingleNode(  std::vector<polyItem> &, std::vector<polyItem> &, std::vector<crossItem> &, Point &, unsigned int, unsigned int, double, double, NextLineTypes);
void pushBackPDoubleNode( std::vector<polyItem> &, std::vector<polyItem> &, std::vector<crossItem> &, Point &, unsigned int, unsigned int, double, NextLineTypes);
void pushBackQDoubleNode( std::vector<polyItem> &, std::vector<polyItem> &, std::vector<crossItem> &, Point &, unsigned int, unsigned int, double, NextLineTypes);
void pushBackPQDoubleNode(std::vector<polyItem> &, std::vector<polyItem> &, std::vector<crossItem> &, Point &, unsigned int, unsigned int, NextLineTypes);
void alterToPDoubleNode( std::vector<polyItem> &, std::vector<polyItem> &, std::vector<crossItem> &, unsigned int, Point &, unsigned int, unsigned int, double, NextLineTypes);
void alterToQDoubleNode( std::vector<polyItem> &, std::vector<polyItem> &, std::vector<crossItem> &, unsigned int, Point &, unsigned int, unsigned int, double, NextLineTypes);
void alterToPQDoubleNode(std::vector<polyItem> &, std::vector<polyItem> &, std::vector<crossItem> &, unsigned int, Point &, unsigned int, unsigned int, NextLineTypes);
void recomputeIndices(unsigned int, unsigned int, unsigned int&, unsigned int&, unsigned int&);


inline std::ostream & operator << (std::ostream &out, const polyItem &poly) {
	out << "[" << ((poly.inside)? "i":"o") << "] cross: ";
	for (unsigned int i = 0; i < poly.inCross.size(); i++) {
		out << std::setw(3) << poly.inCross[i];
	}
	return out;
}

inline std::ostream & operator << (std::ostream &out, const crossItem &cross) {
	out << "(" << std::setw(8) << cross.x << "," << std::setw(8) << cross.y << ") ";
	switch (cross.type) {
	case   Single: out << "  Single(i=" << std::setw(2) << cross.i << ",j=" << std::setw(2) << cross.j << "     "                          << ",m=" << cross.m << ",n=" << std::setw(2) << cross.n <<"     "                          << ")  ";  break;
	case  PDouble: out << " PDouble(i=" << std::setw(2) << cross.i << ",j=" << std::setw(2) << cross.j << ",k=" << std::setw(2) << cross.k << ",m=" << cross.m << ",n=" << std::setw(2) << cross.n <<"     "                          << ")  ";  break;
	case  QDouble: out << " QDouble(i=" << std::setw(2) << cross.i << ",j=" << std::setw(2) << cross.j << "     "                          << ",m=" << cross.m << ",n=" << std::setw(2) << cross.n <<",o=" << std::setw(2) << cross.o << ")  ";  break;
	case PQDouble: out << "PQDouble(i=" << std::setw(2) << cross.i << ",j=" << std::setw(2) << cross.j << ",k=" << std::setw(2) << cross.k << ",m=" << cross.m << ",n=" << std::setw(2) << cross.n <<",o=" << std::setw(2) << cross.o << ")  ";  break;
	case   None:   out << "    None  ";  break;
	}
	out << "pNext=" << ((cross.pNextType==PList)?"p[":"c[") << std::setw(2) << cross.pNext << "] ";
	out << "qNext=" << ((cross.qNextType==QList)?"q[":"c[") << std::setw(2) << cross.qNext << "] ";
	out << "p[" << std::setw(2) << cross.i  << "," << std::setw(8) << cross.pPar  << "] ";
	out << getNextLineTypeName(cross.nextLine) << "  ";
	out << "q[" << std::setw(2) << cross.m  << "," << std::setw(8) << cross.qPar  << "]";
	return out;
}

inline std::ostream & operator << (std::ostream &out, const nodeState &node) {
	out << ((node.nodeType==PList)? "P":((node.nodeType==QList)?"Q":"C")) << "[" << std::setw(2) << node.listInd << "]";
	return out;
}

inline std::ostream & operator << (std::ostream &out, const WeilerAthertonState &state) {
	out << "state " << (state.notEnd?"Continue ":"Ending  ");
	return out;
}

}

#endif /* SRC_MESH_PREPROCESSING_WEILERATHERTON_H_ */

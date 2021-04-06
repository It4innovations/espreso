
#include "weileratherton.h"
#include <iomanip>

using namespace espreso;

CrossTypes espreso::decideCrossType(double pt, double qt, double eps) {
	if (       (0+eps < qt) && (qt < 1-eps)) { // qt   inside
		if (       (0+eps < pt) && (pt < 1-eps)) { // pt  inside
			return Single;
		} else if ((pt < 0-eps) || (1+eps < pt)) { // pt outside
			return None;
		} else {                                   // pt boundary
			return PDouble;
		}
	} else if ((qt <  0-eps) || (1+eps < qt)) {// qt  outside
		return None;
	} else {                                   // qt boundary
		if (       (0+eps < pt) && (pt < 1-eps)) { // pt  inside
			return QDouble;
		} else if ((pt < 0-eps) || (1+eps < pt)) { // pt outside
			return None;
		} else {                                   // pt boundary
			return PQDouble;
		}
	}
}

NextLineTypes espreso::decideSingleCrossNextLineType(  LeftOnRight xmLeftToxixj, LeftOnRight xnLeftToxixj){
	if (           xmLeftToxixj ==  Left && xnLeftToxixj == Right) {
		return PLine;
	} else if (    xmLeftToxixj == Right && xnLeftToxixj ==  Left) {
		return QLine;
	} else {
//		std::cout << "Error(decideSingleCrossNextLineType): xmLeftToxixj=" << getLeftOnRightName(xmLeftToxixj) << " xnLeftToxixj="  << getLeftOnRightName(xnLeftToxixj) << std::endl;
		return NLine;
	}
}

NextLineTypes espreso::decidePDoubleCrossNextLineType( LeftOnRight xmLeftToxixj, LeftOnRight xnLeftToxixj, LeftOnRight xkLeftToxixj, LeftOnRight xmLeftToxjxk, LeftOnRight xnLeftToxjxk, bool isNegativeDot2XjXiXnXm) {
	if (           xmLeftToxixj ==  Left && xnLeftToxixj == Right) {
		if (       xkLeftToxixj ==  Left ) {
			if (   xmLeftToxjxk ==  Left ) {
				return PLine;
			} else {
				return NLine;
			}
		} else {// xkLeftToxixj !=  Left
			if (xnLeftToxjxk ==  Left) {
				return QLine;
			} else {
				return PLine;
			}
		}
	} else if (    xmLeftToxixj == Right && xnLeftToxixj ==  Left) {
		if (       xkLeftToxixj ==  Left ) {
			if (   xnLeftToxjxk ==  Left ) {
				return QLine;
			} else {
				return PLine;
			}
		} else {// xkLeftToxixj !=  Left
			return QLine;
		}
	} else if (    xmLeftToxixj ==    On && xnLeftToxixj ==    On) {
		if (isNegativeDot2XjXiXnXm) { // ij mn oposite direction
			if (       xkLeftToxixj == Right ) {
				return PLine;
			} else {
				return NLine;
			}
		} else {                      // ij mn same    direction
			if (       xkLeftToxixj ==  Left ) {
				return PLine;
			} else {
				return QLine;
			}
		}
	} else {
//		std::cout << "Error(decidePDoubleCrossNextLineType): xmLeftToxixj=" << getLeftOnRightName(xmLeftToxixj) << " xnLeftToxixj="  << getLeftOnRightName(xnLeftToxixj) << " isNegativeDot2XjXiXnXm=" << isNegativeDot2XjXiXnXm << std::endl;
		return NLine;
	}
}

NextLineTypes espreso::decideQDoubleCrossNextLineType( LeftOnRight xmLeftToxixj, LeftOnRight xoLeftToxixj, LeftOnRight xoLeftToxmxn, bool isNegativeDot2XjXiXnXm) {
	if (           xmLeftToxixj ==  Left) {
		if (       xoLeftToxixj ==  Left) {
			if (   xoLeftToxmxn ==  Left ){
				return QLine;
			} else {
				return NLine;
			}
		} else {// xoLeftToxixj !=  Left
			return PLine;
		}
	} else if (    xmLeftToxixj == Right) {
		if (       xoLeftToxixj ==  Left) {
			return QLine;
		} else {// xoLeftToxixj == Right || On
			if (   xoLeftToxmxn == Right ){
				return PLine;
			} else {
				return NLine;
			}
		}
	} else {    // xmLeftToxixj ==    On
		if (isNegativeDot2XjXiXnXm) { // ij mn oposite direction
			if (       xoLeftToxixj ==  Left) {
				return QLine;
			} else {
				return NLine;
			}
		} else {                      // ij mn same    direction
			if (       xoLeftToxixj ==  Left) {
				return QLine;
			} else {
				return PLine;
			}
		}
	}
}

NextLineTypes espreso::decidePQDoubleCrossNextLineType(LeftOnRight xkLeftToxixj, LeftOnRight xoLeftToxixj, LeftOnRight xoLeftToxjxk, LeftOnRight xoLeftToxmxn, LeftOnRight xkLeftToxmxn, LeftOnRight xkLeftToxnxo, LeftOnRight xmLeftToxixj, LeftOnRight xmLeftToxjxk, bool isNegativeDot2XjXiXnXm) {
	if (       xkLeftToxixj == Left) { // angle ijk is in (  0,180) degrees
		if (       xoLeftToxmxn == Left) {                             // angle mno is in (  0,180) degrees
			if (       xoLeftToxixj == Left && xoLeftToxjxk == Left) { //                                 o is  inside angle ijk
				return QLine;
			} else if (xoLeftToxixj == Left && xoLeftToxjxk == On  ) { //                                 o is  on jk
				return PQLine; // PLine
			} else {                                                   //                                 o is outside angle ijk
				if (xkLeftToxmxn == Left  && xkLeftToxnxo == Left) {   //                                            k is  inside angle mno
					return PLine;
				} else {
					return NLine;
				}
			}
		} else if( xoLeftToxmxn == On) {                               // angle mno is         180  degrees
			if (       xoLeftToxixj == Left && xoLeftToxjxk == Left) { //                                 o is  inside angle ijk
				return QLine;
			} else if (xoLeftToxixj == Left && xoLeftToxjxk == On  ) { //                                 o is  on jk
				return PQLine; // PLine
			} else {                                                    //                                o is outside angle ijk
				if (xkLeftToxnxo == Left) {                             //                                           k is  inside angle mno
					return PLine;
				} else {
					return NLine;
				}
			}
		} else{ // xoLeftToxmxn == Right                                  angle mno is in (180,360) degrees
			if (       xoLeftToxixj == Left && xoLeftToxjxk == Left ) { //                                o is  inside angle ijk
				if (   xmLeftToxixj == Left && xmLeftToxjxk == Left ) { //                                           m is  inside angle ijk
					return PmQiLine;
				} else {                                                 //                                          m not inside angle ijk
					return QLine;
				}
			} else if (xoLeftToxixj == Left && xoLeftToxjxk == On   ) { //                                o is  on jk
				return PQLine;
			} else {                                                    //                                o not inside angle ijk
				if (  xmLeftToxixj != Right && xmLeftToxjxk != Right) { //                                           m is  inside angle ijk
					return PLine;
				} else {                                                //                                           m not inside angle ijk
					return NLine;
				}
			}
		}
	} else if( xkLeftToxixj == On ) {  // angle ijk is         180  degrees
		if (       xoLeftToxmxn == Left) {                             // angle mno is in (  0,180) degrees
			if (       xoLeftToxixj == Left ) {                        //                                 o is  inside angle ijk
				return QLine;
			} else if (xoLeftToxixj ==  On  && xmLeftToxixj == Left ) {//                                 o is  on jk
				return PQLine; // PLine
			} else {                                                   //                                 o is outside angle ijk
				if (xkLeftToxmxn == Left  && xkLeftToxnxo == Left) {   //                                            k is  inside angle mno
					return PLine;
				} else {
					return NLine;
				}
			}
		} else if (xoLeftToxmxn == On  ) {                             // angle mno is         180  degrees
			if (       xoLeftToxixj == Left ) {                        //                                 o is  inside angle ijk
				return QLine;
			} else if (xoLeftToxixj ==  On  ) {                        //                                 o is  on jk
				if ( isNegativeDot2XjXiXnXm ) {
					return NLine; // oposite
				} else {
					return PQLine; // PLine // same
				}
			} else {// xoLeftToxixj ==  Right                          //                                 o is outside angle ijk
				return PLine;
			}
		} else{ // xoLeftToxmxn == Right                                  angle mno is in (180,360) degrees
			if (       xoLeftToxixj == Left ) {                        //                                 o is  inside angle ijk
				if (   xmLeftToxixj == Left ) {                        //                                            m is  inside angle ijk
					return PmQiLine;
				} else {                                               //                                            m not inside angle ijk
					return QLine;
				}
			} else if (xoLeftToxixj ==  On  ) {                        //                                 o is  on jk
				if ( isNegativeDot2XjXiXnXm ) {
					return PLine;  // oposite
				} else {
					return PQLine; // same
				}
			} else {// xoLeftToxixj ==  Right                          //                                 o is outside angle ijk
				return PLine;
			}
		}
	} else {// xkLeftToxixj == Right   // angle ijk is in (180,360) degrees
		if (       xoLeftToxmxn == Left) {                             // angle mno is in (  0,180) degrees
			if (       xoLeftToxixj == Left  || xoLeftToxjxk == Left) {//                                 o is  inside angle ijk
				if (       xmLeftToxixj == Left  || xmLeftToxjxk == Left) {//                                        m is  inside angle ijk
					return PmQiLine;
				} else {                                                   //                                        m not inside angle ijk
					return QLine;
				}
			} else if (xoLeftToxixj == Right && xoLeftToxjxk == On  ) {//                                 o is  on jk
				return PQLine;
			} else {                                                   //                                 o is outside angle ijk
				if (   xmLeftToxixj == Left  || xmLeftToxjxk == Left) {//                                            m is  inside angle ijk
					return PLine;
				} else {                                               //                                            m not inside angle ijk
					return NLine;
				}
			}
		} else if (xoLeftToxmxn == On  ) {                             // angle mno is         180  degrees
			if (       xoLeftToxixj == Left  || xoLeftToxjxk == Left) {//                                 o is  inside angle ijk
				if (       xmLeftToxixj == Left  || xmLeftToxjxk == Left) {//                                        m is  inside angle ijk
					return PmQiLine;
				} else {                                                   //                                        m not inside angle ijk
					return QLine;
				}
			} else if (xoLeftToxixj == Right && xoLeftToxjxk == On  ) {//                                 o is  on jk
				return PQLine;
			} else {                                                   //                                 o is outside angle ijk
				return PLine;
			}
		} else {// xoLeftToxmxn == Right                                  angle mno is in (180,360) degrees
			if (       xoLeftToxixj == Left  || xoLeftToxjxk == Left) {//                                 o is  inside angle ijk
				if (       xmLeftToxixj == Left  || xmLeftToxjxk == Left) {//                                        m is  inside angle ijk
					return PmQiLine;
				} else {                                                   //                                        m not inside angle ijk
					return QLine;
				}
			} else if (xoLeftToxixj == Right && xoLeftToxjxk == On  ) {//                                 o is  on jk
				return PQLine;
			} else {                                                   //                                 o is outside angle ijk
				return PLine;
			}
		}
	}
}

void espreso::pushBackSingleNode(std::vector<polyItem> &pList, std::vector<polyItem> &qList, std::vector<crossItem> &cList, Point &intersection, unsigned int i, unsigned int m, double pt, double qt, NextLineTypes nLT) {
	unsigned int j = (i+1) % pList.size();
	unsigned int n = (m+1) % qList.size();
	pList[i].inCross.push_back(cList.size());
	qList[m].inCross.push_back(cList.size());
	cList.push_back(crossItem());
	cList.back().x    = intersection.x;
	cList.back().y    = intersection.y;
	cList.back().i    = i;
	cList.back().j    = j;
	cList.back().m    = m;
	cList.back().n    = n;
	cList.back().pPar = pt;
	cList.back().qPar = qt;
	cList.back().type = Single;
	cList.back().nextLine = nLT;
//	std::cout << "     c[" << std::setw(2) << cList.size()-1 << "] (" << std::setw(3) << cList.back().x << "," << std::setw(3) << cList.back().y << ")[pt,qt]=[" << pt << "," << qt << "], next=" << getNextLineTypeName(nLT) << "   Single added" << std::endl;
}

void espreso::pushBackPDoubleNode(std::vector<polyItem> &pList, std::vector<polyItem> &qList, std::vector<crossItem> &cList, Point &intersection, unsigned int i, unsigned int m, double qt, NextLineTypes nLT) {
	pList[i].inCross.push_back(cList.size());
	qList[m].inCross.push_back(cList.size());
	cList.push_back(crossItem());
	alterToPDoubleNode(pList, qList, cList, cList.size()-1, intersection, i, m, qt, nLT);
//	std::cout << "     c[" << std::setw(2) << cList.size()-1 << "] (" << std::setw(3) << cList.back().x << "," << std::setw(3) << cList.back().y << ")[pt,qt]=[" << 1 << "," << qt << "], next=" << getNextLineTypeName(nLT) << "  PDouble added" << std::endl;
}

void espreso::pushBackQDoubleNode(std::vector<polyItem> &pList, std::vector<polyItem> &qList, std::vector<crossItem> &cList, Point &intersection, unsigned int i, unsigned int m, double pt, NextLineTypes nLT) {
	pList[i].inCross.push_back(cList.size());
	qList[m].inCross.push_back(cList.size());
	cList.push_back(crossItem());
	alterToQDoubleNode(pList, qList, cList, cList.size()-1, intersection, i, m, pt, nLT);
//	std::cout << "     c[" << std::setw(2) << cList.size()-1 << "] (" << std::setw(3) << cList.back().x << "," << std::setw(3) << cList.back().y << ")[pt,qt]=[" << pt << "," << 1 << "], next=" << getNextLineTypeName(nLT) << "  QDouble added" << std::endl;
}

void espreso::pushBackPQDoubleNode(std::vector<polyItem> &pList, std::vector<polyItem> &qList, std::vector<crossItem> &cList, Point &intersection, unsigned int i, unsigned int m, NextLineTypes nLT) {
	pList[i].inCross.push_back(cList.size());
	qList[m].inCross.push_back(cList.size());
	cList.push_back(crossItem());
	alterToPQDoubleNode(pList, qList, cList, cList.size()-1, intersection, i, m, nLT);
//	std::cout << "     c[" << std::setw(2) << cList.size()-1 << "] (" << std::setw(3) << cList.back().x << "," << std::setw(3) << cList.back().y << ")[pt,qt]=[" << 1 << "," << 1 << "], next=" << getNextLineTypeName(nLT) << " PQDouble added" << std::endl;
}

void espreso::alterToPDoubleNode(std::vector<polyItem> &pList, std::vector<polyItem> &qList, std::vector<crossItem> &cList, unsigned int cInd, Point &intersection, unsigned int i, unsigned int m, double qt, NextLineTypes nLT) {
	unsigned int j = (i+1) % pList.size();
	unsigned int k = (i+2) % pList.size();
	unsigned int n = (m+1) % qList.size();
	cList[cInd].x     = intersection.x;
	cList[cInd].y     = intersection.y;
	cList[cInd].i  = i;
	cList[cInd].j  = j;
	cList[cInd].k  = k;
	cList[cInd].m  = m;
	cList[cInd].n  = n;
	cList[cInd].pPar  = 1.0;
	cList[cInd].qPar  = qt;
	cList[cInd].k = k;
	cList[cInd].type  = PDouble;
	cList[cInd].nextLine = nLT;
	pList[j].inside       = true;
	pList[j].insidesetted = true;
}

void espreso::alterToQDoubleNode(std::vector<polyItem> &pList, std::vector<polyItem> &qList, std::vector<crossItem> &cList, unsigned int cInd, Point &intersection, unsigned int i, unsigned int m, double pt, NextLineTypes nLT) {
	unsigned int j = (i+1) % pList.size();
	unsigned int n = (m+1) % qList.size();
	unsigned int o = (m+2) % qList.size();
	cList[cInd].x     = intersection.x;
	cList[cInd].y     = intersection.y;
	cList[cInd].i  = i;
	cList[cInd].j  = j;
	cList[cInd].m  = m;
	cList[cInd].n  = n;
	cList[cInd].o  = o;
	cList[cInd].pPar  = pt;
	cList[cInd].qPar  = 1.0;
	cList[cInd].type  = QDouble;
	cList[cInd].nextLine = nLT;
	qList[n].inside       = true;
	qList[n].insidesetted = true;
}

void espreso::alterToPQDoubleNode(std::vector<polyItem> &pList, std::vector<polyItem> &qList, std::vector<crossItem> &cList, unsigned int cInd, Point &intersection, unsigned int i, unsigned int m, NextLineTypes nLT) {
	unsigned int j = (i+1) % pList.size();
	unsigned int k = (i+2) % pList.size();
	unsigned int n = (m+1) % qList.size();
	unsigned int o = (m+2) % qList.size();
	cList[cInd].x     = intersection.x;
	cList[cInd].y     = intersection.y;
	cList[cInd].i  = i;
	cList[cInd].j  = j;
	cList[cInd].k  = k;
	cList[cInd].m  = m;
	cList[cInd].n  = n;
	cList[cInd].o  = o;
	cList[cInd].pPar  = 1.0;
	cList[cInd].qPar  = 1.0;
	cList[cInd].k = k;
	cList[cInd].o = o;
	cList[cInd].type  = PQDouble;
	cList[cInd].nextLine = nLT;
	pList[j].inside       = true;
	pList[j].insidesetted = true;
	qList[n].inside       = true;
	qList[n].insidesetted = true;
}

void espreso::recomputeIndices(unsigned int j, unsigned int size, unsigned int& i_, unsigned int& j_, unsigned int& k_) {
	j_ = j;
	k_ = (j+1) % size;
	if (j == 0) {
		i_ = size-1;
	} else {
		i_ = (j-1) % size;
	}
}

WeilerAthertonState::WeilerAthertonState(std::vector<polyItem> &p, std::vector<polyItem> &q, std::vector<crossItem> &c):
pList(p), qList(q), cList(c), pVisited(p.size(),false), qVisited(q.size(),false), cVisited(c.size(),false) {
	notEnd            = true;
	processingOutput  = false;
	pVisitedCountdown = p.size();
	cVisitedCountdown = c.size();
	processedNode.nodeType   = PList;
	processedNode.listInd    = 0;
	outputStartNode.nodeType = PList;
	outputStartNode.listInd  = 0;
	previousLineType  = NLine;
}

void WeilerAthertonState::pNext() {
	previousLineType = PLine;
	switch (processedNode.nodeType) {
	case PList:
		if (!pVisited[processedNode.listInd]) {
			pVisited[processedNode.listInd] = true; pVisitedCountdown--;
		}
		if (pList[processedNode.listInd].inCross.size() > 0) {
			processedNode.listInd  = pList[processedNode.listInd].inCross[0];
			processedNode.nodeType = CList;
		} else {
			processedNode.listInd = (processedNode.listInd+1) % pList.size();
		}
		break;
	case QList:
		qVisited[processedNode.listInd] = true;
//		std::cout << "Error pNext() ... cannot be in processedNode.nodeType == QList" << std::endl;
		break;
	case CList:
		if (!cVisited[processedNode.listInd]) {
			cVisited[processedNode.listInd] = true; cVisitedCountdown--;
		}
		switch (cList[processedNode.listInd].type) {
		case   Single:
			break;
		case  QDouble:
			if (!qVisited[cList[processedNode.listInd].n]) { qVisited[cList[processedNode.listInd].n] = true;}
			break;
		case  PDouble:
			if (!pVisited[cList[processedNode.listInd].j]) { pVisited[cList[processedNode.listInd].j] = true; pVisitedCountdown--; }
			break;
		case PQDouble:
			if (!pVisited[cList[processedNode.listInd].j]) { pVisited[cList[processedNode.listInd].j] = true; pVisitedCountdown--; }
			if (!qVisited[cList[processedNode.listInd].n]) { qVisited[cList[processedNode.listInd].n] = true;}
			break;
		default: break;
		}
		processedNode.nodeType = cList[processedNode.listInd].pNextType;
		processedNode.listInd  = cList[processedNode.listInd].pNext;
		break;
	}
//	std::cout << "      pNext()" << std::endl;
}

void WeilerAthertonState::qNext() {
	previousLineType = QLine;
	switch (processedNode.nodeType) {
	case PList:
		if (!pVisited[processedNode.listInd]) { pVisited[processedNode.listInd] = true; pVisitedCountdown--; }
		std::cout << "Error qNext() ... cannot be in processedNode.nodeType == PList" << std::endl;
		break;
	case QList:
		qVisited[processedNode.listInd] = true;
		if (qList[processedNode.listInd].inCross.size() > 0) {
			processedNode.listInd  = qList[processedNode.listInd].inCross[0];
			processedNode.nodeType = CList;
		} else {
			processedNode.listInd = (processedNode.listInd+1) % qList.size();
		}
		break;
	case CList:
		if (!cVisited[processedNode.listInd]) {
			cVisited[processedNode.listInd] = true; cVisitedCountdown--;
		}
		switch (cList[processedNode.listInd].type) {
		case   Single:
			break;
		case  PDouble:
			if (!pVisited[cList[processedNode.listInd].j]) { pVisited[cList[processedNode.listInd].j] = true; pVisitedCountdown--;}
			break;
		case  QDouble:
			if (!qVisited[cList[processedNode.listInd].n]) { qVisited[cList[processedNode.listInd].n] = true;}
			break;
		case PQDouble:
			if (!pVisited[cList[processedNode.listInd].j]) { pVisited[cList[processedNode.listInd].j] = true; pVisitedCountdown--;}
			if (!qVisited[cList[processedNode.listInd].n]) { qVisited[cList[processedNode.listInd].n] = true;}
			break;
		default:
			break;
		}
		processedNode.nodeType = cList[processedNode.listInd].qNextType;
		processedNode.listInd  = cList[processedNode.listInd].qNext;
		break;
	}
//	std::cout << "      qNext()" << std::endl;
}

void WeilerAthertonState::cNext() {
	if (processedNode.nodeType != CList) { std::cout << "Error cNext() called when not in CList node" << std::endl; notEnd = false;}
	bool pGo = false;
	switch (cList[processedNode.listInd].nextLine) {
	case   PQLine:
		switch (cList[processedNode.listInd].pNextType) {
		case PList: if (pList[cList[processedNode.listInd].pNext].inside) {pGo = true;} break;
		case QList: break;
		case CList: pGo = true; break;
		}
		if (pGo) { pNext();
		} else {   qNext(); }
		break;
	case    PLine: pNext(); break;
	case    QLine: qNext(); break;
	case PmQiLine:
		switch (previousLineType) {
		case PLine: qNext(); break;
		default:    pNext(); break;
		}
		break;
	default:
		std::cout << "Error cNext(  NLine)" << std::endl;
		break;
	}
}

void WeilerAthertonState::startProcessOutput() {
	outputStartNode.listInd  = processedNode.listInd;
	outputStartNode.nodeType = processedNode.nodeType;
	processingOutput = true;
	return;
}

bool WeilerAthertonState::canStartProcessOutput() {
	switch (processedNode.nodeType) {
	case PList:
		return true;
		break;
	case QList:
		setEnd();
		std::cout << "Error:  isfalse(state.processingOutput) && state.processedNode.inNodeType == QList" << std::endl;
		return false;
		break;
	default: // CList:
		switch (cList[processedNode.listInd].type) {
		case Single:   return true;  break;
		case PDouble:  return true;  break;
		case QDouble:  return true;  break;
		case PQDouble: return true;  break;
		default:       return false; break;
		}
	}
}

void WeilerAthertonState::endProcessOutput() {
	processingOutput = false;
	return;
}

#ifndef TREECOLLECTION_H
#define TREECOLLECTION_H
#include <vector>
#include <map>
#include <stack>
#include <set>
#include <Eigen/Dense>
#include "utils.hpp"
using namespace boost;

class TreeCollection {

public:
	bool isUniqueSolution(std::vector<std::pair<E, E>> additionalEq);
	std::map<E, double> getSolution(std::map<int, double> coeff, std::vector<std::pair<E, E>> additionalEq) const;
	static TreeCollection constructTree(const graph_t& g, std::set<int>& M, std::set<int>& Mx, std::set<int>& I);

private:
	TreeCollection(int size);
	void addEdge(int source, int target, int dir);
	int getRoot(int vertex) const;
	std::map<E, int> getCharacteristicVector(const E& edge) const;
	std::map<E, double> getPrivateSolution(std::map<int, double> coeff) const;
	static void dfs(const graph_t& g, std::vector<char>& used, std::set<E>& tempEdges, std::set<int>& M, TreeCollection& tree, int v);
	
	std::vector<int> ancestor;
	std::vector<int> direction;
	std::vector<int> depth;
	std::vector<int> tDinast;
	std::set<E> notTreeEdges;
};

#endif // TREECOLLECTION_H

#include "treecollection.hpp"

TreeCollection::TreeCollection(int size) {
	ancestor.resize(size, -1);
	direction.resize(size, 1);
	depth.resize(size);
}

void TreeCollection::addEdge(int source, int target, int dir) {
	ancestor[target] = source;
	depth[target] = depth[source] + 1;
	direction[target] = dir;
}

bool TreeCollection::isUniqueSolution(std::vector<std::pair<E, E>> additionalEq)
{
	std::size_t dim = additionalEq.size();
	std::size_t dimY = notTreeEdges.size();

	if (dimY > dim) {
		return false;
	}
	if (dimY != 0) {
		std::vector<std::map<E, int>> cVectors(dimY);
		int i = 0;
		for (auto e : notTreeEdges) {
			cVectors[i++] = getCharacteristicVector(e);
		}

		Eigen::MatrixXd determinants(dim, dimY);
		for (std::size_t i = 0; i < dim; ++i) {
			for (std::size_t j = 0; j < dimY; ++j) {
				determinants(i, j) = getWithDef(cVectors[j], additionalEq[i].first, 0) - getWithDef(cVectors[j], additionalEq[i].second, 0);
			}
		}
		Eigen::FullPivLU<Eigen::MatrixXd> fpl(determinants);
		return fpl.rank() == dimY;
	}
	else {
		return true;
	}
}

std::map<E, double> TreeCollection::getSolution(std::map<int, double> coeff, std::vector<std::pair<E, E>> additionalEq) const {
	std::map<E, double> solution;
	std::size_t dim = additionalEq.size();
	std::size_t dimY = notTreeEdges.size();
	auto privateSolution = getPrivateSolution(coeff);

	std::vector<std::map<E, int>> cVectors(dimY);
	size_t i = 0;
	for (auto e : notTreeEdges) {
		cVectors[i++] = getCharacteristicVector(e);
	}

	if (dimY != 0) {
		Eigen::VectorXd b(dim);
		Eigen::VectorXd result(dimY);
		Eigen::MatrixXd determinants(dim, dimY);
		for (std::size_t i = 0; i < dim; ++i) {
			for (std::size_t j = 0; j < dimY; ++j) {
				determinants(i, j) = getWithDef(cVectors[j], additionalEq[i].first, 0) -
					getWithDef(cVectors[j], additionalEq[i].second, 0);
			}
			int coeff1 = 1 - notTreeEdges.count(additionalEq[i].first);
			int coeff2 = 1 - notTreeEdges.count(additionalEq[i].second);
			b(i) = -(coeff1*getWithDef(privateSolution, additionalEq[i].first, 0.0) - coeff2*getWithDef(privateSolution, additionalEq[i].second, 0.0));
		}
		result = (determinants.transpose() * determinants).llt().solve(determinants.transpose() * b);
		i = 0;
		for (auto edge : notTreeEdges) {
			solution[edge] = result(i++, 0);
		}
	}
	for (auto it = tDinast.rbegin(); it != tDinast.rend(); ++it) {
		int j = *it;
		E key;
		if (direction[j] == 1) {
			key = std::make_pair(ancestor[j], j);
		}
		else {
			key = std::make_pair(j, ancestor[j]);
		}
		solution[key] = getWithDef(privateSolution, key, 0.0);
		size_t i = 0;
		for (auto edge : notTreeEdges) {
			solution[key] += solution[edge] * getWithDef(cVectors[i++], key, 0);
		}
	}
	return solution;
}

std::map<E, double> TreeCollection::getPrivateSolution(std::map<int, double> coeff) const {
	std::map<E, double> pSolution;
	for (auto it = tDinast.rbegin(); it != tDinast.rend(); ++it) {
		int j = *it;
		double result = -direction[j] * coeff[j];
		if (ancestor[j] != -1) {
			if (direction[j] == 1) {
				pSolution[std::make_pair(ancestor[j], j)] = result;
				coeff[ancestor[j]] -= result;
			}
			else {
				pSolution[std::make_pair(j, ancestor[j])] = result;
				coeff[ancestor[j]] += result;
			}
		}
		else {
			pSolution[std::make_pair(ancestor[j], j)] = result;
		}
	}
	return pSolution;
}

std::map<E, int> TreeCollection::getCharacteristicVector(const E& edge) const {
	std::map<E, int> cVector;
	cVector[edge] = 1;
	int rootI = getRoot(edge.first);
	int rootJ = getRoot(edge.second);
	int i = edge.first;
	int j = edge.second;
	int sign;
	if (rootI == rootJ) {
		int currentDepth;
		if (depth[i] >= depth[j]) {
			currentDepth = depth[j];
			sign = 1;
		}
		else {
			currentDepth = depth[i];
			sign = -1;
			std::swap(i, j);
		}
		while (depth[i] > currentDepth) {
			if (direction[i] == 1) {
				cVector[std::make_pair(ancestor[i], i)] = sign;
			}
			else {
				cVector[std::make_pair(i, ancestor[i])] = -sign;
			}
			i = ancestor[i];
		}
		if (j != edge.second) {
			sign *= -1;
		}
		while (i != j) {
			if (direction[i] == 1) {
				cVector[std::make_pair(ancestor[i], i)] = sign;
			}
			else {
				cVector[std::make_pair(i, ancestor[i])] = -sign;
			}
			if (direction[j] == 1) {
				cVector[std::make_pair(ancestor[j], j)] = -sign;
			}
			else {
				cVector[std::make_pair(j, ancestor[j])] = sign;
			}
			i = ancestor[i];
			j = ancestor[j];
		}
	}
	else {
		while (i != rootI) {
			if (direction[i] == 1) {
				cVector[std::make_pair(ancestor[i], i)] = 1;
			}
			else {
				cVector[std::make_pair(i, ancestor[i])] = -1;
			}
			i = ancestor[i];
		}
		while (j != rootJ) {
			if (direction[j] == 1) {
				cVector[std::make_pair(ancestor[j], j)] = -1;
			}
			else {
				cVector[std::make_pair(j, ancestor[j])] = 1;
			}
			j = ancestor[j];
		}
		cVector[std::make_pair(-1, rootI)] = 1;
		cVector[std::make_pair(-1, rootJ)] = -1;
	}
	return cVector;
}

int TreeCollection::getRoot(int vertex) const {
	while (ancestor[vertex] != -1) {
		vertex = ancestor[vertex];
	}
	return vertex;
}

TreeCollection TreeCollection::constructTree(const graph_t& g, std::set<int>& M, std::set<int>& Mx, std::set<int>& I) {
	std::set<E> tempEdges;
	std::set<std::pair<int, int>> notInTree;
	int n = num_vertices(g);
	TreeCollection tree(n);
	std::vector<char> used(n);

	for (auto v : Mx) {
		auto adjacentVert = adjacent_vertices(v, g);
		for (auto it = adjacentVert.first; it != adjacentVert.second; ++it) {
			tempEdges.insert(std::make_pair(v, *it));
		}
	}

	for (auto v : I) {
		used[v] = true;
	}

	for (auto v : M) {
		used[v] = true;
	}

	for (auto v : I) {
		tree.tDinast.push_back(v);
		dfs(g, used, tempEdges, M, tree, v);
	}

	for (int i = 0; i < used.size(); ++i) {
		if (!used[i]) {
			used[i] = true;
			dfs(g, used, tempEdges, M, tree, i);
		}
	}
	return tree;
}

void TreeCollection::dfs(const graph_t& g, std::vector<char>& used, std::set<E>& tempEdges, std::set<int>& M, TreeCollection& tree, int v) {
	typedef graph_traits < graph_t >::adjacency_iterator adjacency_iterator;
	std::stack<std::pair<std::pair<adjacency_iterator, adjacency_iterator>, int>> dfs_stack;
	dfs_stack.push(std::make_pair(adjacent_vertices(v, g), v));
	while (!dfs_stack.empty()) {
		bool isPush = false;
		auto& topElement = dfs_stack.top();
		int source = topElement.second;
		for (auto& adjacentVert = topElement.first; adjacentVert.first != adjacentVert.second; ++adjacentVert.first) {
			int target = *adjacentVert.first;
			if (!used[target]) {
				if (!tempEdges.count(std::make_pair(source, target))) {
					isPush = true;
					tree.addEdge(source, target, 1);
				}
				else if (!tempEdges.count(std::make_pair(target, source))) {
					isPush = true;
					tree.addEdge(source, target, -1);
				}
				if (isPush) {
					++adjacentVert.first;
					dfs_stack.push(std::make_pair(adjacent_vertices(target, g), target));
					used[target] = true;
					tree.tDinast.push_back(target);
					break;
				}
			}
			else {
				if (!tempEdges.count(std::make_pair(source, target)) &&
					!(tree.direction[source] == -1 && tree.ancestor[source] == target) &&
					!(tree.direction[target] == 1 && tree.ancestor[target] == source)
					) {
					tree.notTreeEdges.insert(std::make_pair(source, target));
				}
			}
		}
		if (!isPush) {
			dfs_stack.pop();
		}
	}
}

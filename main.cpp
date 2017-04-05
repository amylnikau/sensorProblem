#include <fstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <set>
#include <numeric>
#include <boost/range/irange.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/bimap.hpp>
#include <iostream>
#include <omp.h>
#include "treecollection.hpp"
#include "utils.hpp"

using namespace boost;

typedef boost::bimap< int, int > int_bimap;
typedef int_bimap::value_type mapping;

std::set<int> find_sensor_location(const graph_t& g, std::set<int>& I, int lowerBound) {
	int n = num_vertices(g);
	std::vector<int> range(n);
	std::set<int> M;
	bool found = false;

	for (int size = lowerBound; size <= n; ++size) {
		std::iota(range.begin(), range.end(), 0);
		std::vector<int>::iterator r = range.begin() + size;
		do {
			M = std::set<int>(std::begin(range), r);
			std::set<int> Istar;
			std::set_difference(std::begin(I), std::end(I), std::begin(M), std::end(M), std::inserter(Istar, std::end(Istar)));
			auto Mx = getMx(g, M);
			auto IwithoutMstar = getIwithoutMstar(M, Mx, n);
			auto additionalEq = getAdditionalEquations(g, IwithoutMstar);
			int eqCount = n + additionalEq.size() - size;
			int unknownCount = Istar.size();
			for (auto x : IwithoutMstar) {
				unknownCount += out_degree(x, g);
			}
			if (unknownCount > eqCount) {
				continue;
			}

			TreeCollection tree = TreeCollection::constructTree(g, M, Mx, Istar);
			if (tree.isUniqueSolution(additionalEq)) {
				found = true;
				break;
			}

		} while (next_combination(range.begin(), r, range.end()));

		if (found) {
			break;
		}
	}

	return M;

}

std::map<E, double> find_solution(graph_t& g, std::set<int>& I, std::set<int>& M) {
	typename graph_traits<graph_t>::out_edge_iterator ei, ei_end;
	int n = num_vertices(g);
	std::set<int> Istar;
	std::set_difference(std::begin(I), std::end(I), std::begin(M), std::end(M), std::inserter(Istar, std::end(Istar)));
	auto Mx = getMx(g, M);
	auto IwithoutMstar = getIwithoutMstar(M, Mx, n);
	auto additionalEq = getAdditionalEquations(g, IwithoutMstar);

	std::map<E, double> solution;
	for (auto v : M) {
		auto adjacentVert = adjacent_vertices(v, g);
		for (auto it = adjacentVert.first; it != adjacentVert.second; ++it) {
			solution[std::make_pair(v, *it)] = 1;
			solution[std::make_pair(*it, v)] = 1;
		}
		if (I.count(v) != 0) {
			solution[std::make_pair(-1, v)] = 0;
		}
	}
	for (auto v : Mx) {
		for (boost::tie(ei, ei_end) = out_edges(v, g); ei != ei_end; ++ei) {
			auto target = boost::target(*ei, g);
			solution[std::make_pair(v, target)] = 1;
		}
	}

	std::map<int, double> coeff;
	for (auto v : IwithoutMstar) {
		for (boost::tie(ei, ei_end) = out_edges(v, g); ei != ei_end; ++ei) {
			int target = boost::target(*ei, g);
			auto out_key = std::make_pair(v, target);
			coeff[v] -= getWithDef(solution, out_key, 0.0);
			coeff[target] += getWithDef(solution, out_key, 0.0);
		}
	}
	for (auto v : Mx) {
		for (boost::tie(ei, ei_end) = out_edges(v, g); ei != ei_end; ++ei) {
			int target = boost::target(*ei, g);
			auto out_key = std::make_pair(v, target);
			coeff[v] -= getWithDef(solution, out_key, 0.0);
			coeff[target] += getWithDef(solution, out_key, 0.0);
		}
	}
	for (auto v : M) {
		for (boost::tie(ei, ei_end) = out_edges(v, g); ei != ei_end; ++ei) {
			int target = boost::target(*ei, g);
			auto out_key = std::make_pair(v, target);
			coeff[v] -= getWithDef(solution, out_key, 0.0);
			coeff[target] += getWithDef(solution, out_key, 0.0);
		}
	}

	TreeCollection tree = TreeCollection::constructTree(g, M, Mx, Istar);
	auto partSolution = tree.getSolution(coeff, additionalEq);
	solution.insert(partSolution.begin(), partSolution.end());

	return solution;
}

int main() {
	omp_set_num_threads(4);
	int N, edge_num;
	std::set<int> I;
	int_bimap vertexMap;
	std::ifstream input_file("input.txt");
	input_file >> N >> edge_num;
	std::vector<E> edges(edge_num);
	int nodeIndex = 0;
	for (int i = 0; i < edge_num; i += 1)
	{
		int a, b;
		input_file >> a >> b;
		if (vertexMap.left.find(a) == vertexMap.left.end()) {
			vertexMap.insert(mapping(a, nodeIndex++));
		}
		if (vertexMap.left.find(b) == vertexMap.left.end()) {
			vertexMap.insert(mapping(b, nodeIndex++));
		}
		edges[i] = E(vertexMap.left.at(a), vertexMap.left.at(b));
	}
	for (int n; input_file >> n; )
	{
		I.insert(vertexMap.left.at(n));
	}
	input_file.close();

	graph_t g(edges.begin(), edges.end(), N);
	for (int v = 0; v < N; ++v) {
		g[v].vertex_label = vertexMap.right.at(v);
	}
	std::cout << "Our graph" << std::endl;
	print_graph(g, get(&vertex_info::vertex_label, g));

	int lowerBound = find_lower_bound(g, I);
	std::set<int> M = find_sensor_location(g, I, std::max(1, lowerBound));
	std::cout << "\nSensor location using exhaustive search: M = { ";
	for (auto x : M) {
		std::cout << vertexMap.right.at(x) << " ";
	}
	std::cout << "}" << std::endl;
	auto solution = find_solution(g, I, M);
	std::cout << "Solution:" << std::endl;
	std::cout << solution.size() << std::endl;
	for (auto e : solution) {
		if (e.first.first != -1) {
			std::cout << vertexMap.right.at(e.first.first);
		}
		else {
			std::cout << 0;
		}
		std::cout << "," << vertexMap.right.at(e.first.second) << " " << e.second << std::endl;
	}
	return EXIT_SUCCESS;
}

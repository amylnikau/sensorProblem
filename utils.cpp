#include "utils.hpp"

bool operator >(const queue_type& x, const queue_type& y) {
	return std::tie(x.degree, x.node) > std::tie(y.degree, y.node);
}

std::set<int> getMx(const graph_t& g, const std::set<int>& M) {
	/**
		Returns Mplus for a given set M.

		g - Graph.
		M - Set of vertices in which located sensors.
	*/
	graph_traits<graph_t>::adjacency_iterator vi, vi_end;
	std::set<int> Mx;
	std::set<int> tempM;
	for (const auto v : M)
	{
		tie(vi, vi_end) = adjacent_vertices(v, g);
		tempM.insert(vi, vi_end);
	}
	set_difference(tempM.begin(), tempM.end(), M.begin(), M.end(), std::inserter(Mx, Mx.end()));
	return Mx;
}
std::set<int> getIwithoutMstar(const std::set<int>& M, const std::set<int>& Mx, int n) {
	/**
		Returns I\Mstar for a given set M.

		M - Set of vertices in which located sensors.
		Mx - Set of vertices which adjacent to vertices from M.
		n - Number of vertices.
	*/
	std::set<int> Mstar;
	std::set<int> IwithoutMstar;
	std::set<int> I;
	for (int i = 0; i < n; ++i) {
		I.insert(I.end(), i);
	}
	std::set_union(M.begin(), M.end(), Mx.begin(), Mx.end(), std::inserter(Mstar, Mstar.end()));
	set_difference(I.begin(), I.end(), Mstar.begin(), Mstar.end(), std::inserter(IwithoutMstar, IwithoutMstar.end()));
	return IwithoutMstar;
}

std::vector<std::pair<E, E>>
getAdditionalEquations(const graph_t& g, const std::set<int>& IwithoutMstar) {
	/**
		Returns additional equations for a given set M.

		g - Graph.
		IwithoutMstar - I\Mstar.
	*/

	graph_traits<graph_t>::out_edge_iterator ei, ei_end;
	std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> result;
	for (auto x : IwithoutMstar) {
		if (out_degree(x, g) > 1) {
			boost::tie(ei, ei_end) = boost::out_edges(x, g);
			auto baseEdge = std::make_pair(source(*ei, g), target(*ei, g));
			for (++ei; ei != ei_end; ++ei) {
				result.push_back(std::make_pair(std::make_pair(source(*ei, g), target(*ei, g)), baseEdge));
			}
		}
	}
	return result;
}

int find_lower_bound(graph_t& g, std::set<int>& I) {
	typename graph_traits<graph_t>::vertex_iterator vi, vi_end;
	std::vector<int> degrees;
	for (boost::tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi) {
		auto v = vertex(*vi, g);
		degrees.push_back(out_degree(v, g) + I.count(v));
	}
	std::sort(degrees.begin(), degrees.end());
	int sumDegree = 0;
	size_t i = degrees.size() - 1;
	int lowerBound = 0;
	while (sumDegree < I.size()) {
		sumDegree += degrees[i--];
		++lowerBound;
	}
	return lowerBound;
}

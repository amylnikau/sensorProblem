#ifndef UTILS_H
#define UTILS_H
#include <stack>
#include <set>
#include <queue>
#include <vector>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>

using namespace boost;
struct vertex_info {
	int vertex_label;
};
struct queue_type {
	int degree;
	int node;
};

typedef adjacency_list < vecS, vecS, bidirectionalS, vertex_info> graph_t;
typedef std::pair < int, int > E;

std::set<int> getMx(const graph_t& g, const std::set<int>& M);

std::set<int> getIwithoutMstar(const std::set<int>& M, const std::set<int>& Mx, int n);

std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>
getAdditionalEquations(const graph_t& g, const std::set<int>& IwithoutMstar);

int find_lower_bound(graph_t& g, std::set<int>& I);

template <template<class, class, class...> class C, typename K, typename V, typename... Args>
V getWithDef(const C<K, V, Args...>& m, K const& key, const V & defval) {
	typename C<K, V, Args...>::const_iterator it = m.find(key);
	if (it == m.end())
		return defval;
	return it->second;
}

template < class BidirectionalIterator >
bool next_combination(BidirectionalIterator first1,
	BidirectionalIterator last1,
	BidirectionalIterator first2,
	BidirectionalIterator last2) {
	if ((first1 == last1) || (first2 == last2)) {
		return false;
	}
	BidirectionalIterator m1 = last1;
	BidirectionalIterator m2 = last2;
	--m2;
	while (--m1 != first1 && !(*m1 < *m2)) {
	}
	bool result = (m1 == first1) && !(*first1 < *m2);
	if (!result) {
		while (first2 != m2 && !(*m1 < *first2)) {
			++first2;
		}
		first1 = m1;
		std::iter_swap(first1, first2);
		++first1;
		++first2;
	}
	if ((first1 != last1) && (first2 != last2)) {
		m1 = last1;
		m2 = first2;
		while ((m1 != first1) && (m2 != last2)) {
			std::iter_swap(--m1, m2);
			++m2;
		}
		std::reverse(first1, m1);
		std::reverse(first1, last1);
		std::reverse(m2, last2);
		std::reverse(first2, last2);
	}
	return !result;
}

template < class BidirectionalIterator >
bool next_combination(BidirectionalIterator first,
	BidirectionalIterator middle,
	BidirectionalIterator last) {
	return next_combination(first, middle, middle, last);
}

#endif



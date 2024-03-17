/**
 * File: KDTree.h
 * Author: (your name here)
 * ------------------------
 * An interface representing a kd-tree in some number of dimensions. The tree
 * can be constructed from a set of data and then queried for membership and
 * nearest neighbors.
 */

#ifndef KDTREE_INCLUDED
#define KDTREE_INCLUDED

#include "Point.h"
#include "BoundedPQueue.h"
#include <stdexcept>
#include <cmath>

#include <memory>

// "using namespace" in a header file is conventionally frowned upon, but I'm
// including it here so that you may use things like size_t without having to
// type std::size_t every time.
// using namespace std;
using std::size_t;
using std::cout;    using std::endl;
using std::string;  using std::vector;  using std::set;
using std::shared_ptr; using std::make_shared;
using std::nullptr_t;

template <size_t N, typename ElemType>
class KDTree {
public:
    // Constructor: KDTree();
    // Usage: KDTree<3, int> myTree;
    // ----------------------------------------------------
    // Constructs an empty KDTree.
    KDTree();
    
    // Destructor: ~KDTree()
    // Usage: (implicit)
    // ----------------------------------------------------
    // Cleans up all resources used by the KDTree.
    ~KDTree();
    
    // KDTree(const KDTree& rhs);
    // KDTree& operator=(const KDTree& rhs);
    // Usage: KDTree<3, int> one = two;
    // Usage: one = two;
    // -----------------------------------------------------
    // Deep-copies the contents of another KDTree into this one.
    KDTree(const KDTree& rhs);
    KDTree& operator=(const KDTree& rhs);
    
    // size_t dimension() const;
    // Usage: size_t dim = kd.dimension();
    // ----------------------------------------------------
    // Returns the dimension of the points stored in this KDTree.
    size_t dimension() const;
    
    // size_t size() const;
    // bool empty() const;
    // Usage: if (kd.empty())
    // ----------------------------------------------------
    // Returns the number of elements in the kd-tree and whether the tree is
    // empty.
    size_t size() const;
    bool empty() const;
    
    // bool contains(const Point<N>& pt) const;
    // Usage: if (kd.contains(pt))
    // ----------------------------------------------------
    // Returns whether the specified point is contained in the KDTree.
    bool contains(const Point<N>& pt) const;
    
    // void insert(const Point<N>& pt, const ElemType& value);
    // Usage: kd.insert(v, "This value is associated with v.");
    // ----------------------------------------------------
    // Inserts the point pt into the KDTree, associating it with the specified
    // value. If the element already existed in the tree, the new value will
    // overwrite the existing one.
    void insert(const Point<N>& pt, const ElemType& value);
    
    // ElemType& operator[](const Point<N>& pt);
    // Usage: kd[v] = "Some Value";
    // ----------------------------------------------------
    // Returns a reference to the value associated with point pt in the KDTree.
    // If the point does not exist, then it is added to the KDTree using the
    // default value of ElemType as its key.
    ElemType& operator[](const Point<N>& pt);
    
    // ElemType& at(const Point<N>& pt);
    // const ElemType& at(const Point<N>& pt) const;
    // Usage: cout << kd.at(v) << endl;
    // ----------------------------------------------------
    // Returns a reference to the key associated with the point pt. If the point
    // is not in the tree, this function throws an out_of_range exception.
    ElemType& at(const Point<N>& pt);
    const ElemType& at(const Point<N>& pt) const;
    
    // ElemType kNNValue(const Point<N>& key, size_t k) const
    // Usage: cout << kd.kNNValue(v, 3) << endl;
    // ----------------------------------------------------
    // Given a point v and an integer k, finds the k points in the KDTree
    // nearest to v and returns the most common value associated with those
    // points. In the event of a tie, one of the most frequent value will be
    // chosen.
    ElemType kNNValue(const Point<N>& key, size_t k) const;

private:
    // TODO: Add implementation details here.
    struct Node {
        Point<N> key;
        ElemType value;
        size_t level; // level of the node in the KDTree : 0 for x, 1 for y, 2 for z, ...
        shared_ptr<Node> leftc; //! OBS: root & children are *pointers*, as otherwise could not check for nullptr
        shared_ptr<Node> rightc;

        // Add default and parameterized constructor
        Node() : level(0), leftc(nullptr), rightc(nullptr) {};

        Node(Point<N> key, ElemType value, size_t level, nullptr_t leftc, nullptr_t rightc)
        : key(key), value(value), level(level), leftc(leftc), rightc(rightc) {};
    };

    shared_ptr<Node> root; //root of the KDTree 
    size_t numNodes; // number of nodes in the KDTree

    // void deleteNode(Node* node)
    // Usage: deleteNode(root);
    // Recursively deletes the specified node and all its children in the KDTree
    void deleteNode(shared_ptr<Node> node);

    // void rec_ksearch(const Point<N>& key, BoundedPQueue<Point<N>>& bpq, Node* p_current) const
    // Usage: 
    // ----------------------------------------------------
    // Helper function for kNNValue to recurse on tree, filling up the bpq
    void recKSearch(const Point<N>& pt, BoundedPQueue<Point<N>>& bpq, std::shared_ptr<Node> p_current, size_t lvl) const;

    // Node* findNode(const Point<N>& pt) const
    // Usage: Node* foundNode = findNode(v);
    // ----------------------------------------------------
    // Returns a pointer to the node in the KDTree that contains the specified
    // point. If the point is not in the tree, this function returns nullptr.
    shared_ptr<Node> findNode(const Point<N>& pt) const;
  
    // void _copyFromNode(const Node* current)
    // Usage: _copyFromNode(p_node);
    // ----------------------------------------------------
    // Helper function to perform copies of KDTrees, used for both copy constructor and assignment operator
    // Copies the current child_to_copy node and its children, and attaches it to the parent node
    void _copyFromNode(shared_ptr<Node> parent_node, const shared_ptr<Node> child_to_copy, const bool is_lc);    
 };

/** KDTree class implementation details */
/* Constructor : */
template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree() {
    numNodes = 0;
    root = nullptr;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree() {
    deleteNode(root);
    numNodes = 0;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree& rhs) {
    if (this == &rhs) return; //self-assignment
    root = nullptr; //root of the KDTree 
    numNodes = rhs.size();
    if (rhs.empty()) return;
    // copy root & then continue down the tree
    // root = unique_ptr<Node>(new Node {.key = rhs.root->key, .value = rhs.root->value, .level = rhs.root->level, nullptr, nullptr});
    root = make_shared<Node>(rhs.root->key, rhs.root->value,rhs.root->level, nullptr, nullptr);
    _copyFromNode(root, rhs.root->leftc, true);
    _copyFromNode(root, rhs.root->rightc, false);
    return;
}

template <size_t N, typename ElemType>
typename KDTree<N, ElemType>::KDTree& KDTree<N, ElemType>::operator=(const KDTree& rhs){
    if (this == &rhs) return *this; //self-assignment
    deleteNode(root);
    root = nullptr; //root of the KDTree
    numNodes = rhs.size();
    if (rhs.empty()) return *this;
    // copy root & then continue down the tree
    // root = new Node {.key = rhs.root->key, .value = rhs.root->value, .level = rhs.root->level, nullptr, nullptr};
    root = make_shared<Node>(rhs.root->key, rhs.root->value, rhs.root->level, nullptr, nullptr);
    _copyFromNode(root, rhs.root->leftc, true);
    _copyFromNode(root, rhs.root->rightc, false);
    return *this;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const {
    return N;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::size() const {
    return numNodes;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const {
    return numNodes == 0;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N>& pt) const {
    // start at root, and traverse the tree to find the point
    return (findNode(pt)!=nullptr);
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N>& pt, const ElemType& value) {
    // if the tree is empty, create a new node and set it as the root
    if (root == nullptr) {
        // defining the Node the long way: 
        root = shared_ptr<Node>(new Node); //persist this node in the KDTree (past the scope of this function call)
        root->key = pt;
        root->value = value;
        root->level = 0;
        root->leftc = nullptr;
        root->rightc = nullptr;
        numNodes += 1;
        return;
    } else {
        std::shared_ptr<Node> current = root;
        while (current != nullptr) {
            // if the point already exists in the tree, update the value
            if (current->key == pt) {
                current->value = value;
                return;
            }
            size_t i = current->level;
            if (pt[i] > current->key[i]) { // go right
                if (current->rightc != nullptr) {
                    current = current->rightc;
                } else { //insert here and exit
                    current->rightc = make_shared<Node>(pt, value, (i+1)%N, nullptr, nullptr);
                    numNodes += 1;
                    return;
                }
            } else { // go left
                if (current->leftc != nullptr) {
                    current = current->leftc;
                } else { //insert  here and exit
                    //persist this node in the KDTree (past the scope of this function call)                    
                    current->leftc = make_shared<Node>(pt, value, (i+1)%N, nullptr, nullptr);
                    numNodes += 1;
                    return;
                }
            }
        }
    }
}

template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::operator[](const Point<N>& pt) {
    shared_ptr<Node> foundNode = findNode(pt);
    if (foundNode!=nullptr) {
        return foundNode->value;
    } else {
        ElemType default_value = ElemType(); //!
        insert(pt,default_value);
        shared_ptr<Node> insertedNode = findNode(pt);
        return (insertedNode->value);
    }
}

template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::at(const Point<N>& pt) {
    shared_ptr<Node> foundNode = findNode(pt);
    if (foundNode!=nullptr) {
        return foundNode->value;
    } else {
        throw out_of_range("Point was not found.");
    }
}

template <size_t N, typename ElemType>
const ElemType& KDTree<N, ElemType>::at(const Point<N>& pt) const {
    shared_ptr<Node> foundNode = findNode(pt);
    if (foundNode!=nullptr) {
        return foundNode->value;
    } else {
        throw out_of_range("Point was not found.");
    }
}

template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::kNNValue(const Point<N>& key, size_t k) const {
    BoundedPQueue<Point<N>> bpq(k); //! BPQ is of points, not nodes for moment
    recKSearch(key, bpq, root,0);

    // returns the most common value associated with those points, with tie giving most frequent:
    // array of pairs of (Point<N>, ElemType) to store the output of the BPQ
    std::map<ElemType, size_t> freq_map; // value -> frequency
    size_t orig_size = bpq.size();
    for (size_t i = 0; i < orig_size; i++) {
        ElemType pt_value = at(bpq.dequeueMin());
        if (freq_map.count(pt_value)) { //0 or 1
            freq_map[pt_value] += 1;
        } else {
            freq_map[pt_value] = 1;
        }
    }
    ElemType most_freq = std::max_element(freq_map.begin(), freq_map.end(),
        [](const auto& lhs, const auto& rhs) {return lhs.second < rhs.second; })->first;
    return most_freq;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::recKSearch(const Point<N>& pt, BoundedPQueue<Point<N>>& bpq, shared_ptr<Node> p_current, size_t lvl) const {
    //? should the return be bpq instead of void, and bpq be passed-by-reference?
    bool go_left;
    double diff_i;    
    if (p_current == nullptr) {
        return;
    } else {
        // recursively search side of tree the point is on for best fits:
        double priority = Distance(pt,p_current->key);
        bpq.enqueue(p_current->key, priority);
        diff_i = pt[lvl] - p_current->key[lvl];
        go_left = (diff_i < 0) ? true : false;
        // recursively search half of the tree that wud contain test point:
        lvl = (lvl + 1)%N; // update next pt idx to compare
        if (go_left){
            recKSearch(pt, bpq, p_current->leftc, lvl);
        } else {
            recKSearch(pt, bpq, p_current->rightc, lvl);
        }
    }
    // but check that node on other side of splitting hyper-plane isn't closer:
    bool bpq_can_be_filled = (bpq.size() < bpq.maxSize()) && (bpq.size() < numNodes);
    bool crosses_hyperplane = (abs(diff_i) < bpq.worst());
    if (bpq_can_be_filled || crosses_hyperplane) {
        if (go_left) {
            recKSearch(pt, bpq, p_current->rightc, lvl);
        } else {
            recKSearch(pt, bpq, p_current->leftc, lvl); 
        }
    }
    return;        
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::deleteNode(shared_ptr<Node> node) {
        if (node == nullptr) return;
        deleteNode(node->leftc);
        deleteNode(node->rightc);
        // delete node;
        node.reset(); //TODO : shouldn't be necessary
}

template <size_t N, typename ElemType>
shared_ptr<typename KDTree<N, ElemType>::Node> KDTree<N, ElemType>::findNode(const Point<N>& pt) const {
    shared_ptr<Node> p_current = root;
    while (p_current!=nullptr){
        if (p_current->key == pt) {
            return p_current; //Found it
        } else {
            size_t i = p_current->level;
            p_current = (pt[i] > p_current->key[i]) ? p_current->rightc : p_current->leftc;
        }
    }
    return p_current; // return nullptr if not found
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::_copyFromNode(shared_ptr<Node> parent_node, const shared_ptr<Node> child_to_copy, const bool is_lc) {
        if (child_to_copy == nullptr) return;
        shared_ptr<Node> next_node;
        if (parent_node == nullptr) { // if the parent node is nullptr, then the child_to_copy is the root
            root = make_shared<Node>(child_to_copy->key, child_to_copy->value, child_to_copy->level, nullptr, nullptr);
            next_node = root;
        } else { // attach the node to the parent and continue            
            if (is_lc) {
                // parent_node->leftc = new Node {.key = child_to_copy->key, .value = child_to_copy->value, .level = child_to_copy->level, nullptr, nullptr};
                parent_node->leftc = make_shared<Node>(child_to_copy->key, child_to_copy->value, child_to_copy->level, nullptr, nullptr);                
                next_node = parent_node->leftc;
            } else {
                parent_node->rightc = make_shared<Node>(child_to_copy->key, child_to_copy->value, child_to_copy->level, nullptr, nullptr);
                next_node = parent_node->rightc;
            }
        }
        // and continue with the children
        _copyFromNode(next_node, child_to_copy->leftc, true);
        _copyFromNode(next_node, child_to_copy->rightc, false);
    }


#endif // KDTREE_INCLUDED

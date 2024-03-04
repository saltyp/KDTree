/*
Test suite using boosttest.
skip Makefile and simply compile using eg `# g++ --std=c++17 -o test tryBPQ.cpp` for header-only 
or `g++ -std=c++17 tryBPQ.cpp -o egboostcompile /usr/lib/x86_64-linux-gnu/libboost_unit_test_framework.a` for static build usage 
*/

#define BOOST_TEST_MODULE My Test 
#include <boost/test/unit_test.hpp> //using static build
// #include <boost/test/included/unit_test.hpp> //header-only too slow comppilation

#include "Point.h"
#include "BoundedPQueue.h"

using namespace std;

BOOST_AUTO_TEST_CASE(bpq_understand) 
{
    BoundedPQueue<char> bpq(3);
    bpq.enqueue('A', 1);
    BOOST_TEST(bpq.size()==1);
    BOOST_TEST(bpq.maxSize()==3);
    bpq.enqueue('E', 3);
    bpq.enqueue('C', 2);
    bpq.enqueue('Z', 5);
    BOOST_TEST(bpq.size()==3);
    BOOST_TEST(bpq.maxSize()==3);
    BOOST_TEST(!bpq.empty());
    BOOST_TEST(bpq.best()==1,"Best:: highest priority :: lowest number (first in queue): 1");
    BOOST_TEST(bpq.worst()==3,"Worst :: lowest priority :: lowest number (first in queue): 3");
    BOOST_TEST(bpq.dequeueMin()=='A',"We expect dequeueMin to give element with highest priority corresponding to lowest number (first in queue): 1");
    BOOST_TEST(bpq.size()==2,"DequeueMin removes element with lowest number");
    BOOST_TEST(bpq.dequeueMin()=='C',"We expect dequeueMin to give element with next lowest number (next in queue): 2");
    BOOST_TEST(bpq.size()==1,"DequeueMin removes element with lowest number");
    BOOST_TEST(bpq.dequeueMin()=='E',"We expect dequeueMin to give element with next lowest number (next in queue): 3");
    BOOST_TEST(bpq.size()==0,"DequeueMin removes element with lowest number");
}


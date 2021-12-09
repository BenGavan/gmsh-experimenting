//
// Created by Ben Gavan on 12/07/2021.
//

/*
 * To run, need to include flag '-std=c++11' (when running build command (will get an error with vector initializers otherwise))
 */

#include <iostream>
#include <vector>

using namespace std;

//template <function T>
//void printvector(vector<T> v) {
//    for (T x : v) {
//        cout << x << endl;
//    }
//}
//template void printvector<int>(int);

void printvector(vector<int> v) {
    for (int x : v) {
        cout << x << endl;
    }
}

int main() {
    cout << "hey test" << endl;

    vector<int> v = {1,2,3};
    for (int i : v) {
        cout << i << '\n';
    }

    cout << "second vector init test" << endl;

    vector<int> v2 {2,3,4};
    for (int i : v2) {
        cout << i << '\n';
    }

    cout << "third vector test" << endl;
    printvector({9,8,7,6,-1});

    return 0;
}


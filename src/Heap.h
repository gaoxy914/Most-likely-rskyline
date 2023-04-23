#include "Tuple.h"

class Heap {
private:
    vector<pair<double, int> > heap; // ratio, id
    unordered_map<int, int> indexMap; // id -> index
    int heapSize = 0;

    void swap(int i, int j) {
        auto item1 = heap[i];
        auto item2 = heap[j];
        heap[i] = item2;
        heap[j] = item1;
        indexMap[item2.second] = i;
        indexMap[item1.second] = j;
    }

    void adjustup(int child) {
        int parent = (child - 1)/2;
        while (child > 0) {
            if (heap[parent].first < heap[child].first) {
                swap(parent, child);
                child = parent;
                parent = (child - 1)/2;
            } else {
                break;
            }
        }
    }

    void adjustdown(int parent) {
        int child = parent*2 + 1;
        while (child < heapSize) {
            if (child + 1 < heapSize && heap[child].first < heap[child + 1].first) {
                ++ child;
            }
            if (heap[parent].first < heap[child].first) {
                swap(parent, child);
                parent = child;
                child = parent*2 + 1;
            } else {
                break;
            }
        }
    }

    void resign(const pair<double, int>& item) {
        adjustup(indexMap[item.second]);
        adjustdown(indexMap[item.second]);
    }
public:
    Heap() {}
    ~Heap() {}

    void push(const pair<double, int>& item) {
        heap.push_back(item);
        indexMap[item.second] = heapSize;
        adjustup(heapSize++);
    }

    void pop() {
        if (heapSize > 0) {
            auto& item = heap[0];
            swap(0, heapSize - 1);
            heapSize--;
            indexMap.erase(item.second);
            heap.pop_back();
            adjustdown(0);
        }
    }

    void set(const pair<double, int>& item, const pair<double, int>& target) {
        int index = indexMap[item.second];
        heap[index] = target;
        adjustup(index);
        adjustdown(index);
    }

    double get(const int& id) {
        int index = indexMap[id];
        return heap[index].first;
    }

    pair<double, int> top() {
        if (heapSize > 0) {
            return heap[0];
        }
    }

    int size() {
        return heapSize;
    }

    bool empty() {
        return heap.empty();
    }
};
#ifndef OPERATIONS_H
#define OPERATIONS_H

struct OperationStats {
    long long additions = 0;
    long long multiplications = 0;

    void reset() {
        additions = 0;
        multiplications = 0;
    }
};

extern OperationStats stats;

#endif // OPERATIONS_H

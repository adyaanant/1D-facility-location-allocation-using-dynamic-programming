#include <iostream>
#include <vector>
#include <limits>

using namespace std;

const double INF = numeric_limits<double>::max();

// Function to calculate the distance between two points
double distance(double x1, double x2) {
    return abs(x1 - x2);
}

// Dynamic Programming Algorithm
double facilityLocationAllocation(const vector<double>& a, const vector<double>& r, int m, 
                                    vector<int>& optimalLocations) {

    int n = a.size();
    
    // Initialize DP table
    vector<vector<double>> dp(m + 1, vector<double>(n + 1, INF));
    vector<vector<int>> path(m + 1, vector<int>(n + 1, -1));
    
    // Base case
    dp[0][0] = 0;

    // Dynamic programming
    for (int t = 1; t <= m; ++t) {
        for (int i = 1; i <= n; ++i) {
            for (int j = 0; j < i; ++j) {
                double cost = 0;
                for (int k = j + 1; k <= i; ++k) {
                    cost += r[k - 1] * distance(a[k - 1], a[i - 1]);
                }
                if (dp[t][i] > dp[t - 1][j] + cost) {
                    dp[t][i] = dp[t - 1][j] + cost;
                    path[t][i] = j;
                }
            }
        }

        // Print DP table after each iteration
        cout << "DP table after iteration " << t << ":\n";
        for (int row = 0; row <= m; ++row) {
            for (int col = 0; col <= n; ++col) {
                if (dp[row][col] == INF)
                    cout << "INF\t";
                else
                    cout << dp[row][col] << "\t";
            }
            cout << endl;
        }
        cout << endl;
    }

    // Find optimal cost
    double optimalCost = dp[m][n];

    // Backtrack to find optimal solution
    int index = n;
    for (int t = m; t > 0; --t) {
        optimalLocations[t - 1] = index;
        index = path[t][index];
    }

    return optimalCost;
}

int main() {
    vector<double> a = {1, 2, 3, 4, 5, 7, 8, 10}; // Existing facility locations
    vector<double> r = {2.0, 1.0, 2.5, 1.5, 2.5, 4.0, 3.0}; // Requirements at existing facilities
    int m = 3; // Number of variable facilities to be located

    vector<int> optimalLocations(m);
    double optimalCost = facilityLocationAllocation(a, r, m, optimalLocations);

    cout << "Optimal cost: " << optimalCost << endl;
    cout << "Optimal locations for the variable facilities:\n";
    for (int i = 0; i < m; ++i) {
        cout << "Facility " << i + 1 << ": " << a[optimalLocations[i] - 1] << endl;
    }
    return 0;
}


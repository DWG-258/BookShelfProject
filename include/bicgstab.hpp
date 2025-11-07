#ifndef BICGSTAB_H
#define BICGSTAB_H

#include <vector>
#include <iostream>
#include <functional>

class BiCGSTABSolver {
public:
    // 结果结构体
    struct Result {
        std::vector<double> solution;  // 方程解
        int iterations;                // 迭代次数
        double finalResidual;          // 最终残差
        bool converged;                // 是否收敛
        
        void printSolution() const;
        bool verify(const std::vector<std::vector<double>>& A, 
                    const std::vector<double>& b, 
                    double tol = 1e-6) const;
    };
    
    // 构造函数
    BiCGSTABSolver(double tolerance = 1e-8, int maxIterations = 1000);
    
    // 设置参数
    void setTolerance(double tol);
    void setMaxIterations(int maxIter);
    void setVerbosity(bool verbose);
    
    // 求解函数
    Result solve(const std::vector<std::vector<double>>& A, 
                const std::vector<double>& b) const;
    
    // 使用初始猜测求解
    Result solve(const std::vector<std::vector<double>>& A, 
                const std::vector<double>& b,
                const std::vector<double>& initialGuess) const;
                
    // 带回调的求解函数 (用于监控迭代过程)
    Result solve(const std::vector<std::vector<double>>& A, 
                const std::vector<double>& b,
                std::function<void(int, double)> iterationCallback) const;
                
    // 带初始猜测和回调的完整求解函数
    Result solve(const std::vector<std::vector<double>>& A, 
                const std::vector<double>& b,
                const std::vector<double>& initialGuess,
                std::function<void(int, double)> iterationCallback) const;
                
private:
    double m_tolerance;
    int m_maxIterations;
    bool m_verbose;
    
    // 向量运算工具函数
    static double dot(const std::vector<double>& a, const std::vector<double>& b);
    static std::vector<double> matrixVectorProduct(const std::vector<std::vector<double>>& A, 
                                               const std::vector<double>& v);
    static std::vector<double> vectorAdd(const std::vector<double>& a, const std::vector<double>& b);
    static std::vector<double> vectorSub(const std::vector<double>& a, const std::vector<double>& b);
    static std::vector<double> vectorScale(double alpha, const std::vector<double>& a);
    static double vectorNorm(const std::vector<double>& v);
};


#endif // BICGSTAB_H
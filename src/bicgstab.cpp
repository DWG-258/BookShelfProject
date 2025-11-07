#include "bicgstab.hpp"
#include <cmath>
#include "placedata.hpp"

// 构造函数
BiCGSTABSolver::BiCGSTABSolver(double tolerance, int maxIterations)
    : m_tolerance(tolerance), m_maxIterations(maxIterations), m_verbose(true) {}

// 设置参数
void BiCGSTABSolver::setTolerance(double tol) {
    m_tolerance = tol;
}

void BiCGSTABSolver::setMaxIterations(int maxIter) {
    m_maxIterations = maxIter;
}

void BiCGSTABSolver::setVerbosity(bool verbose) {
    m_verbose = verbose;
}

// 向量点乘
double BiCGSTABSolver::dot(const std::vector<double>& a, const std::vector<double>& b) {
    double result = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

// 矩阵乘向量 (A*v)
std::vector<double> BiCGSTABSolver::matrixVectorProduct(
    const std::vector<std::vector<double>>& A, 
    const std::vector<double>& v) {
    
    size_t n = A.size();
    std::vector<double> result(n, 0.0);
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            result[i] += A[i][j] * v[j];
        }
    }
    return result;
}

// 向量加法 (a + b)
std::vector<double> BiCGSTABSolver::vectorAdd(
    const std::vector<double>& a, 
    const std::vector<double>& b) {
    
    size_t n = a.size();
    std::vector<double> result(n);
    
    for (size_t i = 0; i < n; ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

// 向量减法 (a - b)
std::vector<double> BiCGSTABSolver::vectorSub(
    const std::vector<double>& a, 
    const std::vector<double>& b) {
    
    size_t n = a.size();
    std::vector<double> result(n);
    
    for (size_t i = 0; i < n; ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

// 向量数乘 (alpha * a)
std::vector<double> BiCGSTABSolver::vectorScale(
    double alpha, 
    const std::vector<double>& a) {
    
    size_t n = a.size();
    std::vector<double> result(n);
    
    for (size_t i = 0; i < n; ++i) {
        result[i] = alpha * a[i];
    }
    return result;
}

// 计算向量范数
double BiCGSTABSolver::vectorNorm(const std::vector<double>& v) {
    return std::sqrt(dot(v, v));
}

// 基本求解函数
BiCGSTABSolver::Result BiCGSTABSolver::solve(
    const std::vector<std::vector<double>>& A, 
    const std::vector<double>& b) const {
    
    std::vector<double> initialGuess(b.size(), 0.0);
    return solve(A, b, initialGuess);
}

// 带初始猜测的求解函数
BiCGSTABSolver::Result BiCGSTABSolver::solve(
    const std::vector<std::vector<double>>& A, 
    const std::vector<double>& b,
    const std::vector<double>& initialGuess) const {
    
    return solve(A, b, initialGuess, [this](int iter, double residual) {
        if (m_verbose && iter % 10 == 0) {
            std::cout << "迭代 " << iter << ": 残差 = " << residual << std::endl;
        }
    });
}

// 带回调的求解函数
BiCGSTABSolver::Result BiCGSTABSolver::solve(
    const std::vector<std::vector<double>>& A, 
    const std::vector<double>& b,
    std::function<void(int, double)> iterationCallback) const {
    
    std::vector<double> initialGuess(b.size(), 0.0);
    return solve(A, b, initialGuess, iterationCallback);
}

// 带初始猜测和回调的完整求解函数
BiCGSTABSolver::Result BiCGSTABSolver::solve(
    const std::vector<std::vector<double>>& A, 
    const std::vector<double>& b,
    const std::vector<double>& initialGuess,
    std::function<void(int, double)> iterationCallback) const {
    
    size_t n = A.size();
    
    // 检查矩阵维度
    if (n == 0 || A[0].size() != n || b.size() != n) {
        throw std::invalid_argument("矩阵A和向量b的维度不兼容");
    }
    
    // 初始化
    std::vector<double> x = initialGuess;  // 使用初始猜测解
    std::vector<double> r = vectorSub(b, matrixVectorProduct(A, x));  // r0 = b - A*x0
    std::vector<double> r_hat = r;  // 选择r_hat = r0
    
    double rho_prev = 1.0;
    double alpha = 1.0;
    double omega = 1.0;
    
    std::vector<double> v(n, 0.0);
    std::vector<double> p(n, 0.0);
    
    double residualNorm = vectorNorm(r);
    double initialResidual = residualNorm;
    
    // 准备结果
    Result result;
    result.converged = false;
    
    // 迭代主循环
    int iter = 0;
    while (iter < m_maxIterations && residualNorm > m_tolerance * initialResidual) {
        double rho = dot(r_hat, r);
        
        if (std::abs(rho) < 1e-15) {
            if (m_verbose) {
                std::cout << "BiCGSTAB方法失败: rho接近零" << std::endl;
            }
            break;
        }
        
        double beta = (rho / rho_prev) * (alpha / omega);
        
        // p = r + beta * (p - omega * v)
        p = vectorAdd(r, vectorScale(beta, vectorSub(p, vectorScale(omega, v))));
        
        // v = A * p
        v = matrixVectorProduct(A, p);
        
        alpha = rho / dot(r_hat, v);
        
        // s = r - alpha * v
        std::vector<double> s = vectorSub(r, vectorScale(alpha, v));
        
        if (vectorNorm(s) < m_tolerance * initialResidual) {
            // 如果s很小，更新x并退出
            x = vectorAdd(x, vectorScale(alpha, p));
            residualNorm = vectorNorm(s);
            result.converged = true;
            break;
        }
        
        // t = A * s
        std::vector<double> t = matrixVectorProduct(A, s);
        
        omega = dot(t, s) / dot(t, t);
        
        // 更新解向量 x
        x = vectorAdd(vectorAdd(x, vectorScale(alpha, p)), vectorScale(omega, s));
        
        // 更新残差 r
        r = vectorSub(s, vectorScale(omega, t));
        
        rho_prev = rho;
        
        residualNorm = vectorNorm(r);
        ++iter;
        
        // 调用迭代回调
        iterationCallback(iter, residualNorm);
    }
    
    if (iter < m_maxIterations && !result.converged) {
        result.converged = true;
    }
    
    if (m_verbose) {
        if (result.converged) {
            std::cout << "BiCGSTAB在 " << iter << " 次迭代后收敛" << std::endl;
        } else {
            std::cout << "BiCGSTAB未在 " << m_maxIterations << " 次迭代内收敛" << std::endl;
        }
        std::cout << "最终相对残差: " << residualNorm / initialResidual << std::endl;
    }
    
    // 设置结果
    result.solution = x;
    result.iterations = iter;
    result.finalResidual = residualNorm / initialResidual;
    
    return result;
}

// 打印结果
void BiCGSTABSolver::Result::printSolution() const {
    std::cout << "解向量 x = ";
    for (double val : solution) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
    std::cout << "迭代次数: " << iterations << std::endl;
    std::cout << "相对残差: " << finalResidual << std::endl;
    std::cout << "收敛状态: " << (converged ? "已收敛" : "未收敛") << std::endl;
}

// 验证结果
bool BiCGSTABSolver::Result::verify(
    const std::vector<std::vector<double>>& A, 
    const std::vector<double>& b, 
    double tol) const {
    // 计算 A*x
    std::vector<double> Ax = BiCGSTABSolver::matrixVectorProduct(A, solution);
    
    // 计算 ||Ax - b||
    std::vector<double> diff = BiCGSTABSolver::vectorSub(Ax, b);
    double error = BiCGSTABSolver::vectorNorm(diff) / BiCGSTABSolver::vectorNorm(b);
    
    std::cout << "验证结果: ||Ax - b||/||b|| = " << error << std::endl;
    
    std::cout;
}
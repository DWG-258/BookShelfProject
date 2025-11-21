#include "placedata.hpp"

// 对于固定的线长公式，化简为矩阵形式的参数是确定的，所以总线长等于多个线长相加（即不同连线对应的 A, b 累加），最终得到总的线长矩阵参数，最后求导 A, B 参数不变。
// 可以理解为线长公式的参数是固定的，3种情况对总线长的贡献度不同，把每个参数相加，得到总的线长矩阵参数。再求导
// 得到的结果就是A ，b 两个参数矩阵
void caclculate_parm_matrix(
    const vector<shared_ptr<Module>> &modules,
    const vector<shared_ptr<Net>> &nets,
    vector<unordered_map<int, double>> &A,
    vector<double> &b,
    vector<vector<double>>& W_x_off,
     vector<vector<double>>& W_y_off)
{

    int N = modules.size();
    A.assign(N, {});
    b.assign(N, 0.0);
    
    if (W_x_off.size() != N || W_y_off.size() != N){
        throw std::runtime_error("Invalid w_count count");
    }
 

    for (const auto &net : nets)
    {
        const auto &pins = net->netPins;
        int pinCount = pins.size();
        
        if (pinCount < 2)
            continue;


        // 枚举每一对引脚 (p, q)
        for (int i = 0; i < pinCount; ++i)
        {
            for (int j = i + 1; j < pinCount; ++j)
            {
               
               
                auto p = pins[i];
                auto q = pins[j];

                if (p == nullptr || q == nullptr)
                {
                    throw std::runtime_error("Invalid pin pointer");
                }

                shared_ptr<Module> mp = p.get()->module.lock();
                shared_ptr<Module> mq = q.get()->module.lock();

                
                


                if (mp->idx < 0 || mp->idx >= N)
                {
                    throw std::runtime_error("p Invalid module index");
                }
                if (mq->idx < 0 || mq->idx >= N)
                {
                    throw std::runtime_error("q Invalid module index");
                }

                // 考虑偏移（pin相对于module中心的offset）
                double off_p = p->offset.x;
                double off_q = q->offset.x;
                double delta = off_p - off_q;
                double w = W_x_off[mp->idx][mq->idx];

                
                // Case (i): p,q 都可移动
                if (!mp->isFixed && !mq->isFixed)
                {

                    A[mp->idx][mp->idx] += w;
                    A[mq->idx][mq->idx] += w;
                    A[mp->idx][mq->idx] -= w;
                    A[mq->idx][mp->idx] -= w;

                    b[mp->idx] += 2 * w * delta;
                    b[mq->idx] += -2 * w * delta;
                }

                // Case (ii): p movable, q terminal
                if (!mp->isFixed && mq->isFixed)
                {
                    double const_xq = mq->center.x;
                    A[mp->idx][mp->idx] += w;
                    double xq = delta - const_xq;
                    b[mp->idx] += 2 * w * xq;
                }
                // Case (iii): p terminal, q movable
                if (mp->isFixed && !mq->isFixed)
                {
                    double const_xp = mp->center.x;
                    A[mq->idx][mq->idx] += w;
                    double xp = delta + const_xp;
                    b[mq->idx] += 2 * w * xp;
                }
            }
        }
    }
}

std::pair<vector<unordered_map<int, double>>, vector<double>> PlaceData::get_parm_matrix()
{
   
    vector<unordered_map<int, double>> A;
    vector<double> b;
    vector<vector<double>> W_x_off;
    vector<vector<double>> W_y_off;
    computeWeightMatrices(Nets, Nodes, W_x_off, W_y_off);
    caclculate_parm_matrix(Nodes, Nets, A, b,W_x_off,W_y_off);

    return std::make_pair(move(A), move(b));
}


//计算公式
// w_{x,pq}\quad=\frac{1}{P-1}\frac{1}{\left|x_{p}^{\mathrm{pin}}-x_{q}^{\mathrm{pin}}\right|}.\quad
// w_{y,pq}\quad=\frac{1}{P-1}\frac{1}{\left|y_{p}^{\mathrm{pin}}-y_{q}^{\mathrm{pin}}\right|}.\quad
void PlaceData::computeWeightMatrices(
    const vector<shared_ptr<Net>>& nets,
    const vector<shared_ptr<Module>>& modules,
    vector<vector<double>>& W_x_off,
    vector<vector<double>>& W_y_off,
    double eps // 防止除零的小量
){
  int n=modules.size();
  W_x_off.resize(n,vector<double>(n,0.0));
  W_y_off.resize(n,vector<double>(n,0.0));
  for(const auto& net:nets){
    int P=net->netPins.size();
    if(P<2) continue;
    double base = 1.0/(P-1);
    vector<double> Pin_x(P,0.0);
    vector<double> Pin_y(P,0.0);
    for(int p=0;p<P;++p){
      const auto& pin=net->netPins[p];
      const auto& module_ptr=pin->module.lock();
      if(!module_ptr) continue;
      Pin_x[p]=module_ptr->center.x + pin->offset.x;
      Pin_y[p]=module_ptr->center.y + pin->offset.y;
    }

    for(int p=0;p<P;++p){
      const auto& pin_p=net->netPins[p];
      const auto& module_p=pin_p->module.lock();
      if(!module_p) continue;
      int idx_p=module_p->idx;
      for(int q=p+1;q<P;++q){
        const auto& pin_q=net->netPins[q];
        const auto& module_q=pin_q->module.lock();
        if(!module_q) continue;
        int idx_q=module_q->idx; 
        if(idx_p==idx_q) continue;//continue的作用是跳过相同模块
        //相同模块的pin不计算,是因为pin在同一模块上,不会影响布局
        double abs_dx=std::fabs(Pin_x[p]-Pin_x[q]);
        double abs_dy=std::fabs(Pin_y[p]-Pin_y[q]);
        if(abs_dx<eps) abs_dx=eps;
        if(abs_dy<eps) abs_dy=eps;
        double wx=base/abs_dx;
        double wy=base/abs_dy;
        W_x_off[idx_p][idx_q]+=wx;
        W_x_off[idx_q][idx_p]+=wx;
        W_y_off[idx_p][idx_q]+=wy;
        W_y_off[idx_q][idx_p]+=wy;
      }
    }
  }
}
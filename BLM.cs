using System.Collections.Generic;
using Accord.Math;
using Accord.Math.Optimization;

namespace BLMLab
{
    public class BLM
    {
        public BLM()
        {
            
        }
        public BLM(double[,] covMatrix, double[] weight_mkt, double[,] P, double[] Q_views, double riskPremium, double[][,] bounds = null, double[,] omega = null)
        {
            initialize(covMatrix, weight_mkt, P, Q_views, 1 + riskPremium, 1, omega);
            this._bounds = bounds;
        }
        public BLM(double[,] covMatrix, double[] weight_mkt, double[,] P, double[] Q_views, double mktExpectedReturn, double riskFreeRate, double[][,] bounds = null, double[,] omega = null)
        {
            initialize(covMatrix, weight_mkt, P, Q_views, mktExpectedReturn, riskFreeRate, omega);
            this._bounds = bounds;
        }

        private void initialize(double[,] covMatrix, double[] weight_mkt, double[,] P, double[] Q_views, double mktExpectedReturn, double riskFreeRate, double[,] omega = null)
        {
            this._tau = 0.025;
            this._covMatrix = covMatrix; // 资产收益协方差矩阵
            this._weight_mkt = weight_mkt; // 各资产的市值权重
            this._Q_views = Q_views; // 看法矩阵
            this._P = P; // 涉及主观看法的资产矩阵
            this._riskPremium = mktExpectedReturn - riskFreeRate;// 风险溢价
            this.n = _weight_mkt.Length;// 资产数

            // 设置确信度矩阵，如不提供可自行生成（He and Litterman, 1999）
            if (omega == null)
            {
                var L = new List<double>();
                for (var i = 0; i < _P.GetLength(0); i++)
                {
                    var p_i = _P.GetRow(i);
                    var w_k = p_i.Dot(_covMatrix).Dot(p_i.Transpose()).Multiply(_tau);
                    L.Add((double)w_k.GetValue(0));
                }
                this._Omega = Matrix.Diagonal(L.ToArray());
            }
            else
                this._Omega = omega;
        }

        public double[,] Omega
        {
            get { return _Omega; }
            set { _Omega = value; }
        }

        // 计算先验收益：隐含市场均衡收益
        private double[] impliedEquReturn()
        {
            // lambda    
            var _sigma = (double)(_weight_mkt.Dot(_covMatrix).Dot(_weight_mkt.Transpose())).GetValue(0);
            var lda = _riskPremium / _sigma;
            this.lambda = lda;
            // Implied Market Equlibrium Return Vector
            return _covMatrix.Dot(_weight_mkt).Multiply(lda);
        }

        // 计算后验收益：主观观点结合的预期收益
        private double[] newCombinedReturn()
        {
            this.phai = impliedEquReturn();
            var m1 = _covMatrix.Multiply(_tau).Inverse().Add(_P.Transpose().Dot(_Omega.Inverse()).Dot(_P));
            var m2 = _covMatrix.Multiply(_tau).Inverse().Dot(phai).Add(_P.Transpose().Dot(_Omega.Inverse()).Dot(_Q_views));
            return m1.Inverse().Dot(m2);
        }

        // 生成目标函数,及Gradient Function
        private QuadraticObjectiveFunction obj_fun()
        {
            // max w^T * miu - 0.5 * lambda * w^T * Sigma * w
            this.miu = newCombinedReturn();
            var Q = _covMatrix.Multiply(-lambda); // gradient function matrix of quadratic terms
            var d = miu; // vector of linear terms
            return new QuadraticObjectiveFunction(Q, d);
        }

        // 生成约束
        private double[,] constraintMatrix()
        {
            var eqVec = Vector.Create(n, 1.0); //[1,1,...,1] 1*n Vector

            double[][,] matrices;
            if (this._bounds == null)
            {
                var ineqMat = Matrix.Diagonal(n, 1.0); // diag(1,...,1) n*n Matrix
                matrices = new double[][,]
                {
                    Matrix.RowVector(eqVec),
                    ineqMat            
                };//(n+1)*n constraints matrix
            }
            else
            {
                var ineqMatL = Matrix.Diagonal(n, 1.0);// diag(1,...,1) n*n Matrix
                var ineqMatU = Matrix.Diagonal(n, -1.0);// diag(-1,...,-1) n*n Matrix
                matrices = new double[][,]
                {
                    Matrix.RowVector(eqVec),
                    ineqMatL,
                    ineqMatU
                };//(2n+1)*n constraints matrix
            }
            return Matrix.Stack(matrices);
        }

        // 生成约束值
        private double[] constraintValue()
        {
            if (this._bounds == null)
                return Matrix.Concatenate(new double[] { 1.0 }, Vector.Create(n, 0.0));//[1,0,...,0] 1*(n+1) Vector
            else//_bounds:[[b_1,B_1],[b_2,B_2],...,[b_n,B_n]]
            {
                var arr = Matrix.Stack(_bounds);
                var re =  Matrix.Concatenate(new double[][]
                {
                    new double[] { 1.0 },
                    arr.GetColumn(0),
                    arr.GetColumn(1).Multiply(-1)
                });//[1, b_1, ..., b_n, -B_1, ..., -B_n] :1*(2n+1) Vector
                return re;
            }
        }

        // 求解
        public void fit()
        {
            var solver = new GoldfarbIdnani(obj_fun(), constraintMatrix(), constraintValue(), 1);
            // And attempt solve for the max:
            this.success = solver.Maximize();
            this.solution = solver.Solution;
            this.maxValue = solver.Value;
        }



        private double[,] _covMatrix; // 资产收益协方差矩阵
        private double[] _weight_mkt; // 各资产的市值权重
        private double[] _Q_views; // 看法矩阵
        private double[,] _P; // 涉及主观看法的资产矩阵
        private double[,] _Omega;// 确信度矩阵        
        private double _riskPremium; // 风险溢价 = 市场期望回报率 - 无风险利率(mktExpectedReturn - _riskFreeRate)
        private double[][,] _bounds; // 权重约束边界
        
        public double _tau;
        public int n; // 资产的数量
        public double lambda;
        public double[] phai;// 先验收益：隐含市场均衡收益
        public double[] miu;// 后验收益：主观观点结合的预期收益

        public bool success;
        public double[] solution;
        public double maxValue;
    }
}

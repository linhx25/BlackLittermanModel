/// <summary>
/// 资产配置问题中，常常出现协方差矩阵非正定（Non Positive Definete）的问题，由于协方差矩阵是
/// 半正定矩阵，因此原因是其特征值太小，十分接近0。最大可能是由于资产数远大于时间期数，在这种
/// 情况下增加时间期数可尝试解决。如增加样本数后仍无法解决，尝试本解决方法
/// </summary>
namespace RiskParityLab
{
    using Accord.Math;
    using Accord.Statistics;
    using Accord.Math.Decompositions;
    
    public class NonPositiveDefineteMatrixProcessor
    {
        public NonPositiveDefineteMatrixProcessor()
        {

        }
        public NonPositiveDefineteMatrixProcessor(double[,] matrix)
        {
            this._matrix = matrix;
            this.dim = _matrix.GetLength(0);
        }
        public bool isPositiveDefinete(double[,] m)
        {
            return m.IsPositiveDefinite();
        }
        /// <summary>
        /// Higham, Nicholas J. "Computing a nearest symmetric positive semidefinite matrix."Linear algebra and its applications103 (1988): 103-118.
        /// 对非正定矩阵，可以使用此方法获得最接近的对称【半正定】矩阵:
        /// 1.SVD decomposition: E = U * S * V^T
        /// 2.H = V * S * V^T
        /// 3.result = (E + E^T + H + H^T)*0.25
        /// </summary>
        public double[,] NearestSymmetricPositiveSemidefiniteMatrix
        {
            get 
            {
                var m = _matrix;
                var svd = new SingularValueDecomposition(m);
                var (s, v) = (Matrix.Diagonal(svd.Diagonal), svd.RightSingularVectors);
                var H = v.Dot(s).Dot(v.Transpose());
                return (m.Add(m.Transpose()).Add(H).Add(H.Transpose())).Multiply(0.25);
            }            
        }
        /// <summary>
        /// 构造一个对角的协方差矩阵。对角元素为各维度协方差的方差，各维度之间独立
        /// </summary>
        public double[,] DiagonalCovarianceMatrix
        {
            get 
            {                
                var m = Vector.Zeros(dim);
                for (var i = 0; i < dim; i++)
                {
                    var Var = _matrix.GetRow(i).Variance();
                    m[i] = Var;
                }
                return Matrix.Diagonal(m);
            }
        }
        /// <summary>
        /// 给协方差矩阵对角元加上较小的数值，使之正定
        /// </summary>
        public double[,] ModifiedCoviranceMatrix
        {
            get 
            {
                return _matrix.Add(Matrix.Diagonal(Vector.Create(dim, 1e-4)));
            }
        }
        #region Test
        public static void UnitTest()
        {
            var path = @".\TickData\CovMatrix\20160406covMatrix.csv";
            var dt_cov = Test1.read_csv(path);
            var cov = dt_cov.ToMatrix();
            cov = new NonPositiveDefineteMatrixProcessor(cov).ModifiedCoviranceMatrix;
            var pd = cov.IsPositiveDefinite();
            var c = 1;
        }
        #endregion
        private double[,] _matrix;
        private int dim;
    }
}

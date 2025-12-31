package linearalgebra

import (
	"errors"
	"fmt"
	"math"
	"strconv"
)

var (
	MatrixUnlikeDimensionsError = errors.New("Matrices have unlike dimensions")
	OutOfBoundsError            = errors.New("Index is out of bounds")
	IsNotSquareError            = errors.New("Matrix is not square")
	SingularError               = errors.New("Matrix is singular")
)

// Matrix with M rows and N columns
type Matrix struct {
	in [][]float64
}

func (mat *Matrix) AppendColumn(v Vector) error {
	if mat.M() != v.Dim() {
		return VectorUnlikeDimensionsError
	}
	for i, _ := range mat.in {
		val, _ := v.Get(i)
		mat.in[i] = append(mat.in[i], val)
	}
	return nil
}

func (mat *Matrix) AppendRow(v Vector) error {
	if mat.N() != v.Dim() {
		return VectorUnlikeDimensionsError
	}
	mat.in = append(mat.in, v.ToSlice())
	return nil
}

func (mat Matrix) ToReducedRowEchcelon() error {
	for i := 0; i < mat.M(); i++ {
		curr, _ := mat.GetNthRow(i)
		for j := 0; j < i; j++ {
			v, _ := mat.GetNthRow(j)
			val, _ := curr.Get(j)
			v = v.Scale(-1 * val)
			summ, _ := curr.Add(v)
			curr = summ
		}
		piv, _ := curr.Get(i)
		k := i + 1
		for piv == 0 {
			if k >= mat.M() {
				break
			}
			toAdd, _ := mat.GetNthRow(k)
			summ, _ := curr.Add(toAdd)
			curr = summ
			piv, _ = curr.Get(i)
			k++
		}
		now, _ := curr.Get(i)
		if now == 0 {
			mat.SetRow(i, curr)
			return SingularError
		}
		now, _ = curr.Get(i)
		curr = curr.Scale(1 / now)
		mat.SetRow(i, curr)
	}
	return nil
}

func (a Matrix) Copy() Matrix {
	var in [][]float64
	for _, r := range a.in {
		var row []float64
		for _, c := range r {
			row = append(row, c)
		}
		in = append(in, row)
	}
	return MatrixFrom(in)
}

func (a Matrix) BackSub() (Vector, error) {
	mat := a.Copy()
	if err := mat.ToReducedRowEchcelon(); err != nil {
		return Vector{}, err
	}
	var sol []float64
	for i := mat.M() - 1; i >= 0; i-- {
		found := mat.M() - 1 - i
		curr, _ := mat.GetNthRow(i)
		rhs, _ := curr.Get(mat.N() - 1)
		var lhs float64 = 0
		z := 0
		for j := mat.N() - 2; j > mat.N()-2-found; j-- {
			val, _ := curr.Get(j)
			lhs += sol[z] * val
			z++
		}
		sol = append(sol, rhs-lhs)
	}
	var rev []float64
	for i := len(sol) - 1; i >= 0; i-- {
		rev = append(rev, sol[i])
	}
	return VectorFrom(rev), nil
}

func (a Matrix) SetRow(n int, row Vector) error {
	if row.Dim() != a.N() {
		return VectorUnlikeDimensionsError
	}
	if n < 0 || n >= a.M() {
		return OutOfBoundsError
	}
	a.in[n] = row.ToSlice()
	return nil
}

func (a Matrix) SetCol(n int, col Vector) error {
	if col.Dim() != a.M() {
		return VectorUnlikeDimensionsError
	}
	if n < 0 || n >= a.N() {
		return OutOfBoundsError
	}
	for i := 0; i < a.M(); i++ {
		v, _ := col.Get(i)
		a.in[i][n] = v
	}
	return nil
}

func (a Matrix) Inverse() (Matrix, error) {
	d, err := a.Det()
	if err != nil {
		return Matrix{}, err
	}
	if d == 0 {
		return Matrix{}, SingularError
	}
	b, err := a.Adj()
	if err != nil {
		return Matrix{}, err
	}
	return b.ScalarMultiply(1 / d), nil
}

func (a Matrix) Set(i int, j int, val float64) error {
	if i < 0 || i >= a.M() || j < 0 || j >= a.N() {
		return OutOfBoundsError
	}
	a.in[i][j] = val
	return nil
}

func (a Matrix) ScalarMultiply(c float64) Matrix {
	var buf [][]float64
	for _, r := range a.in {
		var x []float64
		for _, v := range r {
			x = append(x, v*c)
		}
		buf = append(buf, x)
	}
	return MatrixFrom(buf)
}

func (a Matrix) Minor(i, j int) (float64, error) {
	if i < 0 || j < 0 || i >= a.M() || j >= a.N() {
		return 0, OutOfBoundsError
	}
	var buf [][]float64
	for x, r := range a.in {
		if x == i {
			continue
		}
		var z []float64
		for y, v := range r {
			if y != j {
				z = append(z, v)
			}
		}
		buf = append(buf, z)
	}
	return MatrixFrom(buf).Det()
}

func (a Matrix) Cofactor(i, j int) (float64, error) {
	d, err := a.Minor(i, j)
	return math.Pow(-1, float64(i+j)) * d, err
}

func (a Matrix) Adj() (Matrix, error) {
	var buf [][]float64
	for i := 0; i < a.M(); i++ {
		var x []float64
		for j := 0; j < a.N(); j++ {
			c, err := a.Cofactor(i, j)
			if err != nil {
				return Matrix{}, nil
			}
			x = append(x, c)
		}
		buf = append(buf, x)
	}
	return MatrixFrom(buf).Transpose(), nil
}

func MatrixFrom(in [][]float64) Matrix {
	return Matrix{in}
}

func (a Matrix) String() string {
	var s string
	var buf []any
	for _, r := range a.in {
		for _, v := range r {
			s += "%9s "
			buf = append(buf, strconv.FormatFloat(v, 'f', 2, 64))
		}
		s += "\n"
	}

	return fmt.Sprintf(s, buf...)
}

func (a Matrix) Det() (float64, error) {
	if !a.IsSquare() {
		return 0, IsNotSquareError
	}
	if a.M() == 1 {
		return a.in[0][0], nil
	}
	if a.M() == 2 {
		return (a.in[0][0] * a.in[1][1]) - (a.in[0][1] * a.in[1][0]), nil
	}
	var s float64
	for i := 0; i < a.N(); i++ {
		var buf [][]float64
		for j, r := range a.in {
			var x []float64
			if j == 0 {
				continue
			}
			for k, v := range r {
				if k != i {
					x = append(x, v)
				}
			}
			buf = append(buf, x)
		}
		d, _ := Matrix{in: buf}.Det()
		s += math.Pow(-1, float64(i)) * d * a.in[0][i]
	}
	return s, nil
}

// A x B
func (a Matrix) Multiply(b Matrix) (Matrix, error) {
	if a.N() != b.M() {
		return Matrix{}, MatrixUnlikeDimensionsError
	}
	var buf [][]float64
	for i := 0; i < a.M(); i++ {
		var x []float64
		for j := 0; j < b.N(); j++ {
			r, _ := a.GetNthRow(i)
			c, _ := b.GetNthCol(j)
			n, _ := r.Dot(c)
			x = append(x, n)
		}
		buf = append(buf, x)
	}
	return MatrixFrom(buf), nil
}

func (a Matrix) GetNthRow(n int) (Vector, error) {
	var buf []float64
	if n >= a.M() || n < 0 {
		return Vector{}, OutOfBoundsError
	}
	buf = a.in[n]
	return VectorFrom(buf), nil
}

func (a Matrix) GetNthCol(n int) (Vector, error) {
	var buf []float64
	if n >= a.N() || n < 0 {
		return Vector{}, OutOfBoundsError
	}
	for i := 0; i < a.M(); i++ {
		buf = append(buf, a.in[i][n])
	}
	return VectorFrom(buf), nil
}

func (a Matrix) Get(i, j int) (float64, error) {
	if i < 0 || i >= a.M() || j < 0 || j >= a.N() {
		return 0, OutOfBoundsError
	}
	return a.in[i][j], nil
}

func IdentityFunc(i int, j int) float64 {
	if i == j {
		return 1
	}
	return 0
}

func CollectColVectors(c []Vector) (Matrix, error) {
	var buf [][]float64
	for i := 0; i < len(c[0].in); i++ {
		var x []float64
		for _, v := range c {
			n, err := v.Get(i)
			if err != nil {
				return Matrix{}, err
			}
			x = append(x, n)
		}
		buf = append(buf, x)
	}
	return MatrixFrom(buf), nil
}

func (a Matrix) Transpose() Matrix {
	var vecs []Vector
	for i := 0; i < a.M(); i++ {
		r, _ := a.GetNthRow(i)
		vecs = append(vecs, r)
	}
	m, _ := CollectColVectors(vecs)
	return m
}

func Identity(n int) Matrix {
	return CreateMatrix(n, n, IdentityFunc)
}

func CreateMatrix(m int, n int, f func(i, j int) float64) Matrix {
	var buf [][]float64
	for x := 0; x < m; x++ {
		var z []float64
		for y := 0; y < n; y++ {
			z = append(z, f(x, y))
		}
		buf = append(buf, z)
	}
	return MatrixFrom(buf)
}

func (a Matrix) IsSameDimension(b Matrix) bool {
	return (a.M() == b.M()) && (a.N() == b.N())
}

func (a Matrix) IsSquare() bool {
	return a.M() == a.N()
}

// Get number or rows of matrix
func (m Matrix) M() int {
	return len(m.in)
}

// Get number of columns of matrix
func (m Matrix) N() int {
	if len(m.in) == 0 {
		return 0
	}
	return len(m.in[0])
}

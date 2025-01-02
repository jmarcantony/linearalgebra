package linearalgebra

import (
	"errors"
	"fmt"
	"math"
	"strings"
)

var (
	VectorUnlikeDimensionsError = errors.New("Vectors have unlike dimensions")
)

type Vector struct {
	in []float64
}

func (v Vector) Dim() int {
	return len(v.in)
}

func VectorFrom(b []float64) Vector {
	return Vector{b}
}

func (v Vector) Get(n int) (float64, error) {
	if n < 0 || n >= len(v.in) {
		return 0, OutOfBoundsError
	}
	return v.in[n], nil
}

func (v Vector) ToMatrix() Matrix {
	var buf [][]float64
	for i := 0; i < v.Dim(); i++ {
		buf = append(buf, []float64{v.in[i]})
	}
	return MatrixFrom(buf)
}

func (v Vector) Mag() float64 {
	x, _ := v.Dot(v)
	return math.Sqrt(x)
}

func (v Vector) Dot(w Vector) (float64, error) {
	var s float64
	if len(v.in) != len(w.in) {
		return 0, VectorUnlikeDimensionsError
	}
	for i := 0; i < len(v.in); i++ {
		s += v.in[i] * w.in[i]
	}
	return s, nil
}

func (v Vector) String() string {
	var x []any
	for _, y := range v.in {
		x = append(x, y)
	}
	return "<" + fmt.Sprintf(strings.Join((strings.Split(strings.Repeat("%f|", v.Dim()), "|"))[:v.Dim()], ", "), x...) + ">"
}

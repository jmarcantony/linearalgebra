package linearalgebra

import (
	"errors"
	"fmt"
)

var (
	UnlikeVariablesError = errors.New("No of varialbes and dimension of vector do not match")
)

type Equation struct {
	c []float64 // Co-efficients of variables
	r float64   // Result on rhs
}

func EquationFrom(v []float64, r float64) Equation {
	return Equation{c: v, r: r}
}

func (e Equation) String() string {
	var buf string
	var a []any
	for i, v := range e.c {
		if i != len(e.c)-1 {
			buf += "%f" + fmt.Sprintf("x%d + ", i+1)
		} else {
			buf += "%f" + fmt.Sprintf("x%d", i+1)
		}
		a = append(a, v)
	}
	buf += fmt.Sprintf(" = %f", e.r)
	return fmt.Sprintf(buf, a...)
}

func (e Equation) Dim() int {
	return len(e.c)
}

func (e Equation) Get(n int) (float64, error) {
	if n < 0 || n >= e.Dim() {
		return 0, OutOfBoundsError
	}
	return e.c[n], nil
}

func (e Equation) Res() float64 {
	return e.r
}

func (e Equation) Plug(v Vector) (float64, error) {
	var s float64
	if e.Dim() != v.Dim() {
		return 0, UnlikeVariablesError
	}
	for i := 0; i < e.Dim(); i++ {
		x, _ := v.Get(i)
		y, _ := e.Get(i)
		s += y * x
	}
	return s, nil
}

func (e Equation) Check(v Vector) (bool, error) {
	s, err := e.Plug(v)
	if err != nil {
		return false, err
	}
	return s == e.r, nil
}

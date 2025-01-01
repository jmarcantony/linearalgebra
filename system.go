package linearalgebra

type System struct {
	e []Equation
}

func (s System) Satisfies(v Vector) (bool, error) {
	for _, e := range s.e {
		if x, err := e.Check(v); err != nil {
			return x, err
		} else if x == false {
			return x, err
		}
	}
	return true, nil
}

func (s System) Dim() int {
	if len(s.e) == 0 {
		return 0
	}
	return s.e[0].Dim()
}

func (s System) Coeff() Matrix {
	var vecs []Vector
	for i := 0; i < s.Dim(); i++ {
		var buf []float64
		for j := 0; j < s.Dim(); j++ {
			e, _ := s.Get(j)
			v, _ := e.Get(i)
			buf = append(buf, v)
		}
		vecs = append(vecs, VectorFrom(buf))
	}
	c, _ := CollectColVectors(vecs)
	return c
}

func (s System) Get(n int) (Equation, error) {
	if n < 0 || n >= s.Dim() {
		return Equation{}, OutOfBoundsError
	}
	return s.e[n], nil
}

func (s System) Res() Vector {
	var buf []float64
	for _, e := range s.e {
		buf = append(buf, e.r)
	}
	return VectorFrom(buf)
}

func SystemFrom(e []Equation) System {
	return System{e: e}
}

func (s System) Solve() (Vector, error) {
	c := s.Coeff()
	i, err := c.Inverse()
	if err != nil {
		return Vector{}, err
	}
	d, err := i.Multiply(s.Res().ToMatrix().Transpose())
	if err != nil {
		return Vector{}, err
	}
	f, _ := d.GetNthCol(0)
	return f, nil
}

package main

import (
	"fmt"
	"linear_algebra/linear_algebra"
	"log"
)

func main() {
	a := linear_algebra.EquationFrom([]float64{2, -3, 4}, 6)
	b := linear_algebra.EquationFrom([]float64{3, 4, -5}, 7)
	c := linear_algebra.EquationFrom([]float64{4, -5, 6}, 8)
	e := linear_algebra.SystemFrom([]linear_algebra.Equation{a, b, c})
	q, err := e.Solve()
	if err != nil {
		log.Fatal(err)
	}
	fmt.Println(q)
	// p, err := e.Satisfies(q)
	// if err != nil {
	// 	log.Fatal(err)
	// }
	f, _ := e.Get(2)
	fmt.Println(f.Check(q))
	fmt.Println(f.Plug(q))
	fmt.Println(f.Res())
	// fmt.Println(p)
}

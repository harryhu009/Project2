//
//  main.cpp
//  project 2
//
//  Created by hu on 2018/11/16.
//  Copyright © 2018年 hu. All rights reserved.
//
#include <iostream>
#include <cmath>
#include "Polynomial.h"
#ifndef MARMOSET_TESTING
int main();
#endif
//start
// Polynomial initialization
void init_poly(poly_t &p, double const init_coeffs[], unsigned int const init_degree)
{
    
    
    if(p.a_coeffs != nullptr)
    {
        delete[] p.a_coeffs;
        p.a_coeffs = nullptr;
        p.degree = 0;
    }
    p.degree = init_degree;
    std::size_t capacity = init_degree + 1;
    p.a_coeffs = new double[capacity];
    for (std::size_t k = 0; k <capacity; k++)
    {
        p.a_coeffs[k] = init_coeffs[k];
    }
    
}
void destroy_poly(poly_t &p)
{
    delete[]p.a_coeffs;
    p.a_coeffs = nullptr;
    p.degree = 0;
}
//Polynomial degree
unsigned int poly_degree(poly_t const &p)
{
    if (p.a_coeffs == nullptr)
    {
        throw 0;
    }
    else
    {
        return p.degree;
    }
}
//Polynomial coefficient
double poly_coeff(poly_t const &p, unsigned int n)
{
    if (p.a_coeffs == nullptr)
    {
        throw 0;
    }
    else
    {
        if (n <= p.degree)
        {
            return p.a_coeffs[n];
        }
        else
        {
            return 0.0;
        }
    }
}
//Polynomial evaluation
double poly_val(poly_t const &p, double const x)
{
    
    if (p.a_coeffs == nullptr)
    {
        throw 0.0;
    }
    else
    {
        double temp = 0.0;
        std::size_t capacity=p.degree+1;
        for (std::size_t k=0;k<=capacity;k++)
        {
            temp+=p.a_coeffs[k]*pow(x, k);
        }
        return temp;
    }
    
}
//Polynomial addition
void poly_add(poly_t &p, poly_t const &q) {
    if (p.a_coeffs == nullptr || q.a_coeffs == nullptr)
    {
        throw 0;
    }
    else
    {
        poly_t temp{nullptr, 0};
        if (p.degree > q.degree)
        {
            temp.degree = p.degree;
            temp.a_coeffs = new double[temp.degree + 1];
            //initialized
            for (std::size_t k = 0; k < temp.degree + 1; k++)
            {
                temp.a_coeffs[k] = 0;
            }
            //pass the value from p array to temp array
            for (std::size_t k = 0; k < q.degree + 1; k++)
            {
                temp.a_coeffs[k] = q.a_coeffs[k];
            }
            //addition
            for (std::size_t k = 0; k < p.degree + 1; k++) {
                p.a_coeffs[k] = temp.a_coeffs[k] + p.a_coeffs[k];
            }
            delete[]temp.a_coeffs;
        }
        else if (p.degree < q.degree)
        {
            temp.degree=q.degree;
            temp.a_coeffs=new double[temp.degree + 1];
            //initialized
            for (std::size_t k = 0; k < temp.degree + 1; k++)
            {
                temp.a_coeffs[k] = 0;
            }
            //pass the value from p array to temp array
            for (std::size_t k = 0; k < p.degree + 1; k++)
            {
                temp.a_coeffs[k] = p.a_coeffs[k];
            }
            destroy_poly(p);
            p.degree = temp.degree;
            std::size_t capacity = p.degree + 1;
            p.a_coeffs = new double[capacity];
            //addition
            for (std::size_t k = 0; k < capacity; k++)
            {
                p.a_coeffs[k] = temp.a_coeffs[k] + q.a_coeffs[k];
            }
            delete[]temp.a_coeffs;
        }
        else {
            temp.degree = p.degree;
            while (p.a_coeffs[temp.degree] + q.a_coeffs[temp.degree] == 0)
            {
                temp.degree--;
            }
            std::size_t capacity = temp.degree + 1;
            //initialized
            temp.a_coeffs = new double[capacity]{0};
            for (std::size_t k = 0; k <capacity; k++)
            {
                temp.a_coeffs[k] = p.a_coeffs[k]+q.a_coeffs[k];
            }
            destroy_poly(p);
            p.degree = temp.degree;
            capacity = p.degree + 1;
            p.a_coeffs=new double[capacity];
            //transfer the value to p.array
            for (std::size_t j = 0; j < capacity; j++)
            {
                p.a_coeffs[j]=temp.a_coeffs[j];
            }
            destroy_poly(temp);
        }
    }
}
//Polynomial subtraction
void poly_subtract(poly_t &p, poly_t const &q)
{
    
    if(p.a_coeffs == nullptr || q.a_coeffs == nullptr)
    {
        throw 0;
    }
    else {
        poly_t temp{nullptr,0};
        if (p.degree > q.degree) {
            temp.degree = p.degree;
            temp.a_coeffs = new double [temp.degree + 1];
            //initialized
            for (std::size_t k = 0; k < temp.degree + 1; k++)
            {
                temp.a_coeffs[k] = 0;
            }
            //pass the value from p array to temp array
            for (std::size_t k = 0; k < q.degree + 1; k++) {
                
                temp.a_coeffs[k] = q.a_coeffs[k];
            }
            //p
            for (std::size_t k = 0; k < p.degree + 1; k++){
                p.a_coeffs[k] = p.a_coeffs[k] - temp.a_coeffs[k];
            }
            delete[]temp.a_coeffs;
        }
        // p.degree < q.degree
        else if (p.degree < q.degree)
        {
            temp.degree = q.degree;
            temp.a_coeffs = new double [temp.degree+ 1];
            //temp initialized
            for (std::size_t k = 0; k < temp.degree+1; k++)
            {
                temp.a_coeffs[k] = 0;
            }
            //pass the value from p array to temp array
            for (std::size_t k = 0; k < p.degree + 1; k++) {
                
                temp.a_coeffs[k] = p.a_coeffs[k];
            }
            destroy_poly(p);
            p.degree=temp.degree;
            std::size_t capacity = p.degree+1;
            //initialized
            p.a_coeffs=new double[capacity]{0};
            //
            for (std::size_t k = 0; k < capacity; k++)
            {
                p.a_coeffs[k] = temp.a_coeffs[k] - q.a_coeffs[k];
            }
            delete[]temp.a_coeffs;
        }
        //if p.degree = q.degree
        else
        {
            temp.degree = p.degree;
            while(p.a_coeffs[temp.degree]-q.a_coeffs[temp.degree]==0 )
            {
                temp.degree--;
            }
            std::size_t capacity = temp.degree+1;
            //initialized for temp to store the substraction
            temp.a_coeffs = new double[capacity]{0};
            for(std::size_t j = 0; j < capacity; j++)
            {
                temp.a_coeffs[j]=p.a_coeffs[j]-q.a_coeffs[j];
            }
            destroy_poly(p);
            //transfer the new degree(temp.degree) to p.degree(new)
            p.degree = temp.degree;
            //initialized for p to store the result for testing
            p.a_coeffs = new double[capacity]{0};
            for(std::size_t j=0; j<capacity; j++)
            {
                p.a_coeffs[j] = temp.a_coeffs[j];
            }
            destroy_poly(temp);
        }
    }
}
//Polynomial multiplication
void poly_multiply(poly_t &p, poly_t const &q)
{
    if(p.a_coeffs == nullptr || q.a_coeffs == nullptr)
    {
        throw 0;
    }
    else
    {
        //when either polynomial equals to zero, the multiplication result is just one zero.
        if((p.a_coeffs[0]==0&&p.degree==0)||(q.a_coeffs[0]==0&&q.degree==0))
        {
            //delete and null p
            destroy_poly(p);
            p.degree=0;
            //let p_array be a new array with one entry for 0
            p.a_coeffs=new double[1]{0};
        }
        else
        {
            poly_t temp{nullptr,0};
            temp.degree = p.degree + q.degree;
            //initialized
            temp.a_coeffs = new double [temp.degree+1]{0};
            for (std::size_t k = 0; k <= q.degree; k++){
                for(std::size_t j = 0; j <= p.degree; j++)
                {
                    temp.a_coeffs[k+j] += p.a_coeffs[j]*q.a_coeffs[k];
                }
            }
            destroy_poly(p);
            p.degree = temp.degree;
            p.a_coeffs = new double [p.degree+1]{0};
            for(std::size_t i = 0; i <= p.degree ; i++)
            {
                p.a_coeffs[i] = temp.a_coeffs[i];
            }
            destroy_poly(temp);
        }
    }
}
//Polynomial division
double poly_divide(poly_t &p, double r)
{
    if (p.a_coeffs == nullptr) {
        throw 0;
    }
    if(p.degree > 0){
        poly_t temp{nullptr,0};
        temp.degree = p.degree;
        temp.a_coeffs=new double[temp.degree + 1];
        int tmp1 = p.degree;
        for (int k = tmp1; k >= 0; --k){
            if(k == tmp1){
                temp.a_coeffs[k] = p.a_coeffs[k];
            }else{
                temp.a_coeffs[k] = p.a_coeffs[k] + temp.a_coeffs[k+1] * r;
            }
        }
        double tmp2=0.0;
        tmp2 = temp.a_coeffs[0];
        destroy_poly(p);
        p.degree=temp.degree - 1;
        p.a_coeffs=new double[p.degree + 1];
        for (int k = p.degree; k >= 0; --k) {
            p.a_coeffs[k] = temp.a_coeffs[k+1];
        }
        destroy_poly(temp);
        return tmp2;
    }else{
        double tmp2=p.a_coeffs[0];
        p.a_coeffs[0]=0;
        return tmp2;
    }
    
}
//Polynomial differentiation
void poly_diff(poly_t &p)
{
    if (p.a_coeffs == nullptr)
    {
        throw 0;
    }
    else
    {
        if (p.degree == 0)
        {
            p.a_coeffs[0]=0;
        }
        else
        {
            poly_t temp{nullptr,0};
            temp.degree = p.degree-1;
            //initialized
            temp.a_coeffs = new double[temp.degree+1]{0};
            //put the value in p to temp
            for (std::size_t k = 0; k < temp.degree+1; k++)
            {
                temp.a_coeffs[k] = p.a_coeffs[k+1]*(k+1);
            }
            destroy_poly(p);
            p.degree=temp.degree;
            //initialized p_array
            p.a_coeffs=new double [p.degree+1]{0};
            for (std::size_t k=0;k<p.degree+1; k++)
            {
                p.a_coeffs[k] = temp.a_coeffs[k];
            }
            for(int i=0;i<p.degree+1;i++)
            {
                std::cout <<p.a_coeffs[i]<<" ";
            }
            destroy_poly(temp);
        }
    }
}
//Polynomial integral approximation
double poly_approx_int(poly_t const &p, double a, double b, unsigned int n)
{
    if(p.a_coeffs == nullptr)
    {
        throw 0;
    }
    else
    {
        poly_t temp{nullptr,0};
        temp.degree = p.degree;
        //initialized
        temp.a_coeffs = new double[temp.degree+1]{0};
        for(std::size_t k=0;k<temp.degree+1;k++)
        {
            temp.a_coeffs[k]=p.a_coeffs[k];
        }
        double h = (b-a)/n;
        double tmp1 = 0.0;
        double tmp2 = 0.0;
        double tmp3 = 0.0;
        for(std::size_t k=1;k<=n-1;k++)
        {
            double x = a+k*h;
            for (std::size_t j=0;j<=temp.degree;j++)
            {
                tmp1+=temp.a_coeffs[j]*pow(x, j);
            }
        }
        for (std::size_t j=0;j<=temp.degree;j++)
        {
            tmp2+=temp.a_coeffs[j]*pow(a, j);
        }
        for (std::size_t j=0;j<=temp.degree;j++)
        {
            tmp3+=temp.a_coeffs[j]*pow(b, j);
        }
        destroy_poly(temp);
        return (tmp2+tmp3)*h/2+h*tmp1;
    }
}
#ifndef MARMOSET_TESTING
int main() {
    
    poly_t data{nullptr, 0};
    double d[5]{10, 18, 19,33,23};
    init_poly(data, d, 4);
    poly_t data2{nullptr, 0};
    double d2[2]{208.243, 0.180908};
    init_poly(data2, d2, 1);
    
    std::cout<<poly_divide(data2,3.05124);
    std::cout<<std::endl;
    std::cout<<"answer: "<<std::endl;
    
    /*
     poly_diff(data);
     std::cout<<std::endl;
     std::cout<<poly_approx_int(data, -4,4,8) <<std::endl;
     std::cout << "val测试：";
     std::cout << std::endl;
     std::cout << poly_val(data, 3) << std::endl;
     std::cout << "第一个"<<std::endl;
     for (std::size_t k{0}; k < data.degree + 1; ++k) {
     std::cout << data.a_coeffs[k] << " ";
     }
     std::cout<<std::endl;
     std::cout << "第二个"<<std::endl;
     for (std::size_t k{0}; k < data2.degree + 1; ++k) {
     std::cout << data2.a_coeffs[k] << " ";
     }
     std::cout << std::endl;
     std::cout << "第一个*第一个";
     std::cout << std::endl;
     poly_multiply(data, data);
     std::cout << std::endl;
     std::cout << "第一个*第二个";
     std::cout << std::endl;
     poly_multiply(data, data2);
     std::cout << "第一个*第二个";
     std::cout << std::endl;
     */
    return 0;
}
#endif



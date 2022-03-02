#pragma once
#include<iostream>
#include<exception>

class MatrixProblems {
public:
    virtual char const* what() const noexcept = 0;
};

class ZeroMN : public MatrixProblems {
public:
    char const* what() const noexcept override
    {
        return "Error: Matrix cannot have 0 rows or 0 columns!";
    }
};

class NotSquare : public MatrixProblems {
public:
    char const* what() const noexcept override
    {
        return "Error: Matrix is not a square matrix!";
    }
};

class NoInverse : public MatrixProblems {
public:
    char const* what() const noexcept override
    {
        return "Error: Inverse matrix doesn't exist!";
    }
};

class ConcatProblem: public MatrixProblems {
public:
    char const* what() const noexcept override
    {
        return "Error: Matrix should have same amount of rows for concatenation!";
    }
};

class SizeProblem: public MatrixProblems {
public:
    char const* what() const noexcept override
    {
        return "Error: Matrix should have same sizes!";
    }
};

class NotAVector : public MatrixProblems {
public:
    char const* what() const noexcept override
    {
        return "Error: Operand should be a vector!";
    }
};

class RowsAndCols : public MatrixProblems {
public:
    char const* what() const noexcept override
    {
        return "Error: All rows' sizes shoud equal M. All columns' sizes shoud equal N.";
    }
};

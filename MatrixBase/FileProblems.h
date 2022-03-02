#pragma once
#include<iostream>
#include<exception>

class FileProblem {
public:
    virtual char const* what() const noexcept = 0;
};

class OpenProblem : public FileProblem {
public:
    char const* what() const noexcept override
    {
        return "Error: Cannot open the file!";
    }
};

class FormatProblem : public FileProblem {
public:
    char const* what() const noexcept override
    {
        return "Error: Format of file is incorrect!";
    }
};

class UnknownExtention : public FileProblem {
public:
    char const* what() const noexcept override
    {
        return "Error: Extention of file is incorrect!";
    }
};
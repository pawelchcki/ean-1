2+2
using Gadfly, MPFI
A =  Interval[ 1.0 -1 2 -1 -8
       2 -2 3 -3 -20
       1  1 1  0  -2
       1 -1 4  3  4 ]

import Base.isless, Base.promote_rule, Base.convert

promote_rule(::Type{Interval}, ::Type{Int64}) = Interval
promote_rule(::Type{Interval}, ::Type{Float64}) = Interval


# Base.promote_type(Interval, Float64)
# promote_rule(::Type{Int64}, ::Type{Interval}) = Interval
typeof([2])
promote(Interval("1.1"), 1)

Interval("1.1") + Int64[1]

function isless(a::Interval, b::Interval)
  true
end


function reduce_augmented_matrix(A)
  n = size(A,1)
  X = zeros(n)
  X[n] = A[n,n+1]/A[n,n]
  for i = n-1:-1:1
    X[i] = (A[i,n+1] - sum([A[i,j]*X[j] for j = i+1:n ])) / A[i,i]
  end
  return X
end

function subst(A, i, p)
  if p != i
    t = A[i,:]
    A[i,:] = A[p,:]
    A[p,:] = t
  end
end

function reduce_row(A, i)
  n = size(A,1)
  for j = i+1:n
    m = A[j,i]/A[i,i]
    A[j,:] = A[j,:] - A[i,:]*m
  end
end


function gaussian_backwards_substitution(A)
  A=copy(A)
  n = size(A,1)
  for i = 1:n-1
    p = findfirst(A[i:n,i]) # find first nonzero
    if p == 0
      println("no unique solution exists")
      return 0
    end
    p = p + i-1 # normalize coordinates

    subst(A, p, i)
    reduce_row(A,i)

    if A[n,n] == 0
      println("no unique solution exists")
      return 0;
    end
  end

  return reduce_augmented_matrix(A);
end

function gaussian_partial_pivoting(A)
  A = copy(A)
  n = size(A,1)
  for i = 1:n-1
    abs_array = [abs(A[j,i]) for j = i:n]
    p = findfirst(abs_array, maximum(abs_array)) + i-1

    if A[p, i] == 0
      println("no unique solution exists")
      return 0
    end

    subst(A, p, i)
    reduce_row(A, i)
  end

  return reduce_augmented_matrix(A)
end


X = gaussian_partial_pivoting(A)

# println(A[:,1:n] \ A[:,n+1])
println(X)
with_bigfloat_precision(50) do
  println(gaussian_partial_pivoting(A))
end



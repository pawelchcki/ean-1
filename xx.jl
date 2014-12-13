2+2
using Gadfly, MPFI
A =  BigFloat[ 1.0 -1 2 -1 -8
       2 -2 3 -3 -20
       1  1 1  0  -2
       1 -1 4  3  4 ]

function reduce_augmented_matrix(A)
  n = size(A,1)
  X = zeros(eltype(A), n)
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
  Aint = convert(Array{Interval, 2}, A)
  println(Aint);
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

    subst(Aint, p, i)
    reduce_row(Aint, i)

    if A[n,n] == 0
      println("no unique solution exists")
      return 0;
    end
  end
  println("x: ", reduce_augmented_matrix(Aint))
  return reduce_augmented_matrix(A);
end

function gaussian_partial_pivoting(A)
  A = copy(A)
  n = size(A,1)
  for i = 1:n-1
    abs_array = [abs(A[j,i]) for j = i:n]
    p = findfirst(abs_array, maximum(abs_array)) + i-1
    println(abs_array)
    if A[p, i] == 0
      println("no unique solution exists")
      return 0
    end

    subst(A, p, i)
    reduce_row(A, i)
  end

  return reduce_augmented_matrix(A)
end


with_bigfloat_precision(50) do
  println(gaussian_backwards_substitution(A))
end

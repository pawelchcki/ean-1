2+2
using Gadfly
A =  Float64[ 1 -1 2 -1 -8
       2 -2 3 -3 -20
       1  1 1  0  -2
       1 -1 4  3  4 ]
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

    if p != i
      t = A[i,:]
      A[i,:] = A[p,:]
      A[p,:] = t
    end

    for j = i+1:n
      m = A[j,i]/A[i,i]
      A[j,:] = A[j,:] - A[i,:]*m
    end
    if A[n,n] == 0
      println("no unique solution exists")
      return 0;
    end
  end
  X[n] = A[n,n+1]/A[n,n]
  for i = n-1:-1:1
    X[i] = (A[i,n+1] - sum([A[i,j]*X[j] for j = i+1:n ])) / A[i,i]
  end

  return X;
end

gaussian_backwards_substitution(A)

rp = [i for i = 1:n]
println(rp)
for i = 1:n-1
  abs_array = [abs(A[j,i]) for j = i:n]
  println(abs_array)
  val = maximum(abs_array)
  println(val)

  p = findfirst(abs_array, val) + i-1

  if A[p, i] == 0
    println("no unique solution exists")
    return 0
  end

  if p != i
    t = A[i,:]
    A[i,:] = A[p,:]
    A[p,:] = t
  end

  for j = i+1:n
    m = A[j,i] / A[i,i]
    A[j,:] = A[j,:] - m*A[i,:]
  end
end

println(A)


using   MPFI
import Base: convert, isless, promote_rule
function isless(x::Interval, y::Interval)
  return ccall((:mpfi_cmp_default, :libmpfi), Int32, (Ptr{Interval}, Ptr{Interval}), &x, &y) == -1
end
convert(::Type{Interval}, x::ASCIIString) = Interval(x)
convert(::Type{BigFloat}, x::ASCIIString) = BigFloat(x)
promote_rule(::Type{Interval}, ::Type{Int64}) = Interval

A =  [ "1.0" -1.0 2 -1 -8
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
prepareA(A) = convert(Array{BigFloat, 2}, copy(A))
prepareA_interval(A) = convert(Array{Interval, 2}, A)

function gaussian_backwards_substitution(A, interval=false)
  Aint = prepareA_interval(A)
  A = prepareA(A)

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
  if interval
    return reduce_augmented_matrix(Aint)
  else
    return reduce_augmented_matrix(A);
  end
end

function gaussian_partial_pivoting(A)
  A = prepareA_interval(A)

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

set_bigfloat_precision(80)
# A = BigFloat[30.0 591400 591700
#              5.291 -6.130 46.78]


# println(gaussian_partial_pivoting(A))
X = prepareA_interval(A)

println(gaussian_partial_pivoting(A))


A = X[:,1:4]
b = X[:, 5]
import Base.real
real(x::Interval) = x
println(A\b)

#验证所推导的求傅里叶系数公式正确性
#比较variance大小
using LinearAlgebra
using StatsBase
using Plots
using Notifier

function pt(arg...)
    flush(stdout)
    println(arg...)
end
pt()

M=4
m=9
p=normalize!(rand(M),1)#态对应的概率幅list
function A(x)
    sum(p .* cos.(Array(1:M) .* 2x))/2+0.5
end
xk=Array(1:m-1) .* pi ./ m
function est2(i)
    4/m*(1+sum(A.(xk) .* cos.(2i .* xk)))
end
function on(a,b)
    a%b==0 ? b-1 : -1
end
function est3(i)
    #跟est2本质一样，但舍弃est2中的浮点运算
    oo=0
    for j=1:M
        oo+=p[j]*(on(i+j,m)+on(i-j,m))
    end
    4/m*(1+1/2*on(i,m)+1/4*oo)
end
for i=1:M
    pt()
    est2(i)|>pt
    p[i]|>pt
    est3(i)|>pt
end
function D1(N)
    o=0
    for i=1:M
        o+=16/m^2/N*sum(A.(xk) .* (1 .- A.(xk)) .* cos.(2i .* xk).^2)
    end
    o
end
function D2(N)
    o=0
    for i=1:M
        o+=p[i]*(1-p[i])/N/(m-1)
    end
    o
end
pt("variance")
D1(3)|>pt
D2(3)|>pt
# plot(D1.(1:100))
# plot!(D2.(1:100))
# for i=1:m,j=1:m
#     pt(i," ",j)
#     sum(cos.(2i .* Array(1:m-1) .* pi ./ m) .* cos.(2j .* Array(1:m-1) .* pi ./ m))|>pt
# end

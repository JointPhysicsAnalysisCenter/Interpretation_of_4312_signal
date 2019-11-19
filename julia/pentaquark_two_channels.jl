using QuadGK
using Plots

# conversion
const mD0 = 1.86483; const mD0sq = mD0^2
const mJψ = 3.0969; const mJψsq = mJψ^2
const mK = 0.493677; const mKsq = mK^2;
const mp = 0.938; const mpsq = mp^2
const mΛb = 5.6196; const mΛbsq = mΛb^2
const mΣc = 2.4529; const mΣcsq = mΣc^2
const eth0 = mΣc+mD0 ; const sth0 = eth0^2

# conversion between the energy from threshold (in MeV) to s (in GeV)
sofe(e; eth::Real=eth) = (1e-3*e+eth)^2
eofs(s; eth::Real=eth) = (sqrt(s)-eth)*1e3
#
sofe_ΣcD0(e) = sofe(e; eth=eth0)
eofs_ΣcD0(s) = eofs(s; eth=eth0)

# general fumtulas
poly(s,ps) = sum(s^(i-1)*b for (i,b) in enumerate(ps))
mmatrix_amp(s, M, iρ) = inv(M(s) - iρ(s))
# Σc D0
mmatrix_amp_Jψp_ΣcD0_sqrt(s, M) = mmatrix_amp(s, M, s->[-sqrt((mJψ+mp)^2-s) 0.0im
                                                          0.0im             -sqrt((mΣc+mD0)^2-s)])
mmatrix_amp_Jψp_ΣcD0_sqrt_III(s, M) = mmatrix_amp(s, M, s->[sqrt((mJψ+mp)^2-s) 0.0im
                                                            0.0im           -sqrt((mΣc+mD0)^2-s)])
mmatrix_amp_Jψp_ΣcD0_sqrt_IV(s, M) = mmatrix_amp(s, M, s->[-sqrt((mJψ+mp)^2-s) 0.0im
                                                            0.0im            sqrt((mΣc+mD0)^2-s)])
# scattering length approximation
mmatrix_amp_Jψp_ΣcD0_sqrt_scattlen(s, c) = mmatrix_amp_Jψp_ΣcD0_sqrt(s, s->[c[1] c[3]
                                                                            c[3] c[2]])
# scattering amplitude
amp0(s) = mmatrix_amp_Jψp_ΣcD0_sqrt_scattlen(s, (c11=2.6, c22 = 0.22, c12 = 0.85))
amp11(s) = amp0(s)[1,1]
# second sheet amplitude
amp11_III(s) = let c = (c11=2.6, c22 = 0.22, c12 = 0.85)
    mmatrix_amp_Jψp_ΣcD0_sqrt_III(s, s->[c.c11 c.c12
                                        c.c12 c.c22])[1,1]
end
amp11_IV(s) = let c = (c11=2.6, c22 = 0.22, c12 = 0.85)
    mmatrix_amp_Jψp_ΣcD0_sqrt_IV(s, s->[c.c11 c.c12
                                         c.c12 c.c22])[1,1]
end

# production amplitude
prod(s) = poly(s, (p0 = 423.16, p1 = -23.53))
prodamp(s) = prod(s)*amp11(s)

# background
bg(s) = poly(s, (b0 = 402.95, b1 = -15))
λ(x,y,z) = x^2+y^2+z^2-2*x*y-2*y*z-2*z*x
phsp(s) = sqrt(λ(s,mΛbsq,mKsq)*λ(s,mpsq,mJψsq))/(4*sqrt(s))
intensity(s) = (abs2(prodamp(s+1e-6im))+bg(s)) * phsp(s)

# convolution
gaussian(Δe,σ=2.7) = exp(-(Δe)^2/(2*σ^2))/(sqrt(2π)*σ) # in MeV
conv(f; σ=2.7) = e->quadgk(ep->f(ep)*gaussian(e-ep,σ), e-3σ, e+3σ)[1] # in MeV

# LHCb
σe_LHCb(e) = 2.71-6.56e-6*(e-4567)^2 # in MeV
conv_LHCb(f; σ=2.7) = e->quadgk(ep->f(ep)*gaussian(e-ep,σe_LHCb(ep)), e-3σ, e+3σ)[1] # in MeV

# convolution
intens_conv = conv(e->intensity(sofe_ΣcD0(e)))
# check that it works
@show intens_conv(eth0)

pyplot() # the plots are nicer
let ev = range(-50,50,length=200)
    sv = sofe_ΣcD0.(ev)
    #
    calv = prodamp.(sv .+ 1e-6im)
    #
    plot(layout=grid(1,2), size=(950,350), color_palette=palette(:wong2), xlab="e-eth (MeV)")
    plot!(sp=1, ev, [real.(calv) imag.(calv)], frame=:origin, title="amplitude", lab=["re" "im"])
    #
    calv = intens_conv.(ev)
    plot!(sp=2, frame=:origin, title="intensity", ylab="rate (Events / 2 MeV)")
    plot!(sp=2, ev, calv, lab="", lw=1.5)
    calv = map(s->bg(s)*phsp(s), sv)
    plot!(sp=2, ev, calv, lab="", l=(2,0.9), fill_between=fill(0.0,length(calv)), α=0.2,
        ann=(25,440,text("other PWs", 15)))
    calv = map(s->abs2(prodamp(s+1e-6im))*phsp(s), sv)
    plot!(sp=2, ev, calv, lab="", lw=1.5, ann=(25,50,text("studied PW", 15)))
end
savefig(joinpath("plots", "Pc4312.pdf"))

let
    ev = range(-3,3,length=200); sv = sofe_ΣcD0.(ev)
    isv = range(-0.02,0.02,length=100)
    #
    calv = [i>0 ? amp11(r+1im*i) : amp11_III(r+1im*i) for (i,r) in Iterators.product(isv, sv)]
    #
    heatmap(ev, isv, log.(abs2.(calv)), c=:viridis, size=(500,350))
    contour!(ev, isv, log.(abs2.(calv)))
    #
    isv_cut = isv[isv.>0]
    ev_cut = ev[ev.>0]; sv_cut = sofe_ΣcD0.(ev_cut)
    calv = [amp11_IV(r+1im*i) for (i,r) in Iterators.product(isv_cut, sv_cut)]
    contour!(eofs_ΣcD0.(sv_cut), isv_cut, log.(abs2.(calv)), levels=10)
    plot!(xlab = "sqrt(Re(s))-eth (MeV)", ylab = "Im(s) (GeV)")
end
savefig(joinpath("plots", "pole.pdf"))

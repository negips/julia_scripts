

using Revise
using PolynomialBases
using Blink, Plots

# define nodal bases
p = 5 # polynomial degree
basis1 = LobattoLegendre(p)
basis2 = GaussLegendre(p)

# the function that will be interpolated
ufunc(x) = sinpi(x); uprim(x) = π*cospi(x)
#ufunc(x) = 1 / (1 + 25x^2); uprim(x) = -ufunc(x)^2*50x

for basis in (basis1, basis2)
    u = ufunc.(basis.nodes)

    xplot = range(-1, stop=1, length=500)
    uplot = interpolate(xplot, u, basis)

    fig1 = plot(xplot, ufunc.(xplot), label="u", xguide="x", yguide="u")
    plot!(fig1, xplot, uplot, label="I(u)")

    fig2 = plot(xplot, uprim.(xplot), label="u'", xguide="x", yguide="u'")
    plot!(fig2, xplot, interpolate(xplot, basis.D*u, basis), label="I(u)'")

    display(basis)
    display(plot(fig1, fig2))
end


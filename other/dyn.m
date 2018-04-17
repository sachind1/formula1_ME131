function dvdt = dyn(t,v,u,ab)
dvdt = ab(1)*v + ab(2)*u(t);
end

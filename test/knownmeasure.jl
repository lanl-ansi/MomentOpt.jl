@testset Lebesgue measure begin
    @polyvar x y
    m = LebesgueMeasure(variable_box(x => [0,1]))
    @test sprint(show, m) isa String 
    @test integrate.(x.^collect(0:9), ml) ≈ collect(1:10).^(-1)

end

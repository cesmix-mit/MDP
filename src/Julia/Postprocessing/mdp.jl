function mdp(app)

# preprocess input data
app, config = Preprocessing.preprocessing(app); 

# generate code
Gencode.gencode(app); # Generate C++ code

# compile code
Gencode.compile(app);

# run code
Gencode.runcode(app);

#cd(app.currentdir);
if app.training > 0
    filename = app.sourcepath * "C++/Main/coefficients.bin";
    tmp = reinterpret(Float64,read(filename));
    app.coeff = reshape(tmp, 1, length(tmp));
end

return app, config

end

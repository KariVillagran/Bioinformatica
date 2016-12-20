var express 	= require('express');
var app 		= express();
var http 		= require('http');
var server 		= http.createServer(app);
var io 			= require('socket.io').listen(server);
var bodyParser 	= require('body-parser');
var rio 		= require("rio");
var formidable = require('formidable');

var allowCrossDomain = function(req, res, next) {
    res.header('Access-Control-Allow-Origin', "*");
    res.header('Access-Control-Allow-Methods', 'GET,PUT,POST,DELETE');
    res.header('Access-Control-Allow-Headers', 'Content-Type');
    next();
}

app.use(allowCrossDomain); 

// Generamos el puerto de conexión
var port = 8083;

server.listen(port);
console.log('Servidor Node.js iniciado en puerto: ' + port);



// support json encoded bodies
app.use(bodyParser.json()); 
// support encoded bodies
app.use(bodyParser.urlencoded({ extended: true })); 

// Implementacion de metodos de R

app.get('/tree', function (req, res) {
	try{
			var fileName = req.query.fileName;
			var typeFile = req.query.typeFile;
			var model = req.query.model;
			var typeFileUpper = typeFile.toUpperCase();
			var type_method = { 
			  njhamming: `nj(dist.hamming(data))`,
			  EvolutionaryModelHamming: `nj(dist.${typeFile}(as.${typeFileUpper}bin(data),model='${model}'))`,
			  BioNJHamming: `bionj(dist.hamming(data))`,
			  EvolutionaryModelBioNJ: `bionj(dist.${typeFile}(as.${typeFileUpper}bin(data),model='${model}'))`,
			  UPGMAHamming: `upgma(dist.hamming(data))`,
			  UPGMAEvolutionaryModel: `upgma(dist.${typeFile}(as.${typeFileUpper}bin(data),model='${model}'))`,
			  WPGMAHamming: `wpgma(dist.hamming(data))`,
			  WPGMAEvolutionaryModel: `wpgma(dist.${typeFile}(as.${typeFileUpper}bin(data),model='${model}'))`,
			  FastMinimumEvolutionBalanced: `fastme.bal(dist.hamming(data), nni = TRUE, spr = TRUE, tbr = TRUE)`,
			  FastMinimumEvolutionOLS: `fastme.ols(dist.hamming(data), nni = TRUE)`,
			  LeastSquaresFm: `Rfitch(dist.hamming(data), path='C:/Users/Ana/Documents/Usach/Bioinformatica/Proyecto_filogenia/node_services/sources/exe',quiet=TRUE,method='fm')`,
			  LeastSquaresLs: `Rfitch(dist.hamming(data), path='C:/Users/Ana/Documents/Usach/Bioinformatica/Proyecto_filogenia/node_services/sources/exe',quiet=TRUE,method='ls')`
			};
			var type_methode = type_method[req.query.type_m];
		  	rio.$e({
			    command: `library("phangorn");library("ape");library("RJSONIO");library("Rphylip");
			    		  data=as.phyDat(read.${ typeFile }("http://localhost:8081/${ fileName }"));
			    		  write.tree(${type_methode},"C:/Users/Ana/Documents/Usach/Bioinformatica/Proyecto_filogenia/node_services/public/data.tree");
			    		  response <- read.table("C:/Users/Ana/Documents/Usach/Bioinformatica/Proyecto_filogenia/node_services/public/data.tree");
			    		  toString(response$V1)`
			}).then(function (data) {
				res.json(data);
			}).catch(function (err) {
			    console.log(err);
			    res.json(err);
			});
	}
	catch(error){
		throw error;
	}
});


// Arboles de calidad

app.get('/treeBootstrap', function (req, res) {
	try{
		var fileName = req.query.fileName;
		var typeFile = req.query.typeFile;
		var model = req.query.model;
		var type_method = { 
			  njhamming: `nj(dist.hamming(data))`,
			  EvolutionaryModelHamming: `nj(dist.${typeFile}(as.${typeFile}bin(data),model='${model}'))`,
			  BioNJHamming: `bionj(dist.hamming(data))`,
			  EvolutionaryModelBioNJ: `bionj(dist.${typeFile}(as.${typeFile}bin(data),model='${model}'))`,
			  UPGMAHamming: `upgma(dist.hamming(data))`,
			  UPGMAEvolutionaryModel: `upgma(dist.${typeFile}(as.${typeFile}bin(data),model='${model}'))`,
			  WPGMAHamming: `wpgma(dist.hamming(data))`,
			  WPGMAEvolutionaryModel: `wpgma(dist.${typeFile}(as.${typeFile}bin(data),model='${model}'))`,
			  FastMinimumEvolutionBalanced: `fastme.bal(dist.hamming(data), nni = TRUE, spr = TRUE, tbr = TRUE)`,
			  FastMinimumEvolutionOLS: `fastme.ols(dist.hamming(data), nni = TRUE)`,
			  LeastSquaresFm: `Rfitch(dist.hamming(data), path='${__dirname}/sources/exe',quiet=TRUE,method='fm')`,
			  LeastSquaresLs: `Rfitch(dist.hamming(data), path='${__dirname}/sources/exe',quiet=TRUE,method='ls')`
			};
		var type_methode = type_method[req.query.type_m];
	  	rio.$e({
		    command: `library("phangorn");library("ape");library("RJSONIO");library("Rphylip");
			    	  data=as.phyDat(read.${typeFile}("http://localhost:8081/${fileName}"));
			    	  tree1=${type_methode};
  					  NJtrees = bootstrap.phyDat(data, FUN=function(x)${type_methode}, bs=100);
  					  treeNJ <- plotBS(tree1, NJtrees, "phylogram");
  					  write.tree(treeNJ,"C:/Users/Ana/Documents/Usach/Bioinformatica/Proyecto_filogenia/node_services/public/data.tree");
			    	  response <- read.table("C:/Users/Ana/Documents/Usach/Bioinformatica/Proyecto_filogenia/node_services/public/data.tree");
			    	  toString(response$V1)`

		}).then(function (data) {
		    res.json(data);
		}).catch(function (err) {
		    console.log(err);
		    res.json(err);
		});
	}
	catch(error){
		throw error;
	}
});


app.get('/treeParsimony', function (req, res) {
	try{
		var fileName = req.query.fileName;
		var typeFile = req.query.typeFile;
		var model = req.query.model;
		var type_method = { 
			  njhamming: `nj(dist.hamming(data))`,
			  EvolutionaryModelHamming: `nj(dist.${typeFile}(as.${typeFile}bin(data),model='${model}'))`,
			  BioNJHamming: `bionj(dist.hamming(data))`,
			  EvolutionaryModelBioNJ: `bionj(dist.${typeFile}(as.${typeFile}bin(data),model='${model}'))`,
			  UPGMAHamming: `upgma(dist.hamming(data))`,
			  UPGMAEvolutionaryModel: `upgma(dist.${typeFile}(as.${typeFile}bin(data),model='${model}'))`,
			  WPGMAHamming: `wpgma(dist.hamming(data))`,
			  WPGMAEvolutionaryModel: `wpgma(dist.${typeFile}(as.${typeFile}bin(data),model='${model}'))`,
			  FastMinimumEvolutionBalanced: `fastme.bal(dist.hamming(data), nni = TRUE, spr = TRUE, tbr = TRUE)`,
			  FastMinimumEvolutionOLS: `fastme.ols(dist.hamming(data), nni = TRUE)`,
			  LeastSquaresFm: `Rfitch(dist.hamming(data), path='${__dirname}/sources/exe',quiet=TRUE,method='fm')`,
			  LeastSquaresLs: `Rfitch(dist.hamming(data), path='${__dirname}/sources/exe',quiet=TRUE,method='ls')`
			};
		var type_methode = type_method[req.query.type_m];
	  	rio.$e({
		    command: `library("phangorn");library("ape");library("RJSONIO");library("Rphylip");
			    	  data=as.phyDat(read.${typeFile}("http://localhost:8081/${fileName}"));			    	  
			    	  tree1=${type_methode};
			    	  tree1=optim.parsimony(tree1, data);
					  NJtrees = bootstrap.phyDat(data, FUN=function(x)optim.parsimony(tree1,x), bs=100);
					  tree1=acctran(tree1,data);
					  treeNJ <- plotBS(tree1, NJtrees, "phylogram");
					  write.tree(treeNJ,"C:/Users/Ana/Documents/Usach/Bioinformatica/Proyecto_filogenia/node_services/public/data.tree");
			    	  response <- read.table("C:/Users/Ana/Documents/Usach/Bioinformatica/Proyecto_filogenia/node_services/public/data.tree");
			    	  toString(response$V1)`
		}).then(function (data) {
		    res.json(data);
		}).catch(function (err) {
		    console.log(err);
		    res.json(err);
		});
	}
	catch(error){
		throw error;
	}
});

app.get('/treeLikelihood', function (req, res) {
	try{
		var fileName = req.query.fileName;
		var typeFile = req.query.typeFile;
		var model = req.query.model;
		var type_method = { 
			  njhamming: `nj(dist.hamming(data))`,
			  EvolutionaryModelHamming: `nj(dist.${typeFile}(as.${typeFile}bin(data),model='${model}'))`,
			  BioNJHamming: `bionj(dist.hamming(data))`,
			  EvolutionaryModelBioNJ: `bionj(dist.${typeFile}(as.${typeFile}bin(data),model='${model}'))`,
			  UPGMAHamming: `upgma(dist.hamming(data))`,
			  UPGMAEvolutionaryModel: `upgma(dist.${typeFile}(as.${typeFile}bin(data),model='${model}'))`,
			  WPGMAHamming: `wpgma(dist.hamming(data))`,
			  WPGMAEvolutionaryModel: `wpgma(dist.${typeFile}(as.${typeFile}bin(data),model='${model}'))`,
			  FastMinimumEvolutionBalanced: `fastme.bal(dist.hamming(data), nni = TRUE, spr = TRUE, tbr = TRUE)`,
			  FastMinimumEvolutionOLS: `fastme.ols(dist.hamming(data), nni = TRUE)`,
			  LeastSquaresFm: `Rfitch(dist.hamming(data), path='${__dirname}/sources/exe',quiet=TRUE,method='fm')`,
			  LeastSquaresLs: `Rfitch(dist.hamming(data), path='${__dirname}/sources/exe',quiet=TRUE,method='ls')`
			};
		var type_methode = type_method[req.query.type_m];
	  	rio.$e({
		    command: `library("phangorn");library("ape");library("RJSONIO");library("Rphylip");
			    	  data=as.phyDat(read.${typeFile}("http://localhost:8081/${fileName}"));			    	  
			    	  tree1=${type_methode};
			    	  
			    	  mT = NULL;
  					  mT = modelTest(data,tree1);
  					  env <- attr(mT, "env");
  					  ev_tree <- eval(get(mT$Model[which.min(mT$AIC)], env), env);

  					  if ((mT$Model[which.min(mT$AIC)]=="GTR") | (mT$Model[which.min(mT$AIC)]=="GTR+I") | (mT$Model[which.min(mT$AIC)]=="GTR+G")| (mT$Model[which.min(mT$AIC)]=="GTR+G+I")) {ev_tree$model="GTR"};
  					  if ((mT$Model[which.min(mT$AIC)]=="HKY") | (mT$Model[which.min(mT$AIC)]=="HKY+I") | (mT$Model[which.min(mT$AIC)]=="HKY+G")| (mT$Model[which.min(mT$AIC)]=="HKY+G+I")) {ev_tree$model="HKY"};
  					  if ((mT$Model[which.min(mT$AIC)]=="JC") | (mT$Model[which.min(mT$AIC)]=="JC+I") | (mT$Model[which.min(mT$AIC)]=="JC+G")| (mT$Model[which.min(mT$AIC)]=="JC+G+I")) {ev_tree$model="JC"};

  					  tree1=pml(tree = ev_tree$tree, data = data, bf = ev_tree$bf, Q = ev_tree$Q, inv = ev_tree$inv,k = ev_tree$k, shape = ev_tree$shape);

  					  bs <- bootstrap.pml(tree1, bs=100, optNni=TRUE);
  					  treeBS <- plotBS(tree1$tree,bs,"phylogram");

					  write.tree(treeBS,"C:/Users/Ana/Documents/Usach/Bioinformatica/Proyecto_filogenia/node_services/public/data.tree");
			    	  response <- read.table("C:/Users/Ana/Documents/Usach/Bioinformatica/Proyecto_filogenia/node_services/public/data.tree");
			    	  toString(response$V1)`
		}).then(function (data) {
		    res.json(data);
		}).catch(function (err) {
		    console.log(err);
		    res.json(err);
		});
	}
	catch(error){
		throw error;
	}
});



// Subida de archivos al repositorio local

app.post('/upload', function (req, res){
    var form = new formidable.IncomingForm();

    form.parse(req);

    form.on('fileBegin', function (name, file){
        file.path = __dirname + '/public/' + file.name;
    });

    form.on('file', function (name, file){
        console.log('Uploaded ' + file.name);
    });

    res.json('{ok}');
});

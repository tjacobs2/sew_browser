/*
 * Module dependencies
 */
var express = require('express')
  , logger = require('morgan')
  , fs = require('fs');


var app = express();
app.use(express.static(__dirname + '/public')) //set static file location as the public directory
app.use(logger('dev')) //Log every request

// api ---------------------------------------------------------------------
app.get('/graph', function(req, res) {
  var graph_data = fs.readFileSync('sew.json');
  var json_data = JSON.parse(graph_data);
  res.json(json_data);
});

app.get('/models/:model_id', function(req, res) {
  var model_id = req.params.model_id;
  var pdb = fs.readFileSync('pdbs/'+model_id+'.pdb');
  res.send(pdb);
});

// application -------------------------------------------------------------
app.all('/', function(req, res) {
  res.sendfile('./public/index.html');
});

app.listen(3000)

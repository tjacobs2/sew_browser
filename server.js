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
  // console.log(json_data);
  res.json(json_data);
});

// application -------------------------------------------------------------
app.all('/', function(req, res) {
  res.sendfile('./public/index.html');
});

app.listen(3000)

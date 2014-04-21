var align_and_combine = function(protein1, resnum1, protein2, resnum2) {
  var Pcx, Pcy, Pcz, Qcx, Qcy, Qcz;

  var stub1 = protein1.get_residue_stub(resnum1);
  var stub2 = protein2.get_residue_stub(resnum2);

  var coords1 = protein1.get_coord_matrix();
  var coords2 = protein2.get_coord_matrix();

  //Find the center of mass for stub1 and stub2, translate coordinates
  Pcx = numeric.sum(_.map(stub1, function(coord){ return coord[0]; }))/stub1.length;
  Pcy = numeric.sum(_.map(stub1, function(coord){ return coord[1]; }))/stub1.length;
  Pcz = numeric.sum(_.map(stub1, function(coord){ return coord[2]; }))/stub1.length;
  stub1 = _.map(stub1, function(coord){ return [ coord[0]-Pcx, coord[1]-Pcy, coord[2]-Pcz ] });

  Qcx = numeric.sum(_.map(stub2, function(coord){ return coord[0]; }))/stub2.length;
  Qcy = numeric.sum(_.map(stub2, function(coord){ return coord[1]; }))/stub2.length;
  Qcz = numeric.sum(_.map(stub2, function(coord){ return coord[2]; }))/stub2.length;
  stub2 = _.map(stub2, function(coord){ return [ coord[0]-Qcx, coord[1]-Qcy, coord[2]-Qcz ] });

  //Find rotation matrix from stub2 -> stub1 
  //var rot_matrix = kabsch(stub1, stub2);

  var svd_res = numeric.svd(numeric.dot(numeric.transpose(stub2), stub1));
  var V = svd_res.U;
  var S = svd_res.S;
  var Wt = numeric.transpose(svd_res.V);
  var d = (numeric.det(V) * numeric.det(Wt)) < 0.0
  console.log(d);
  if(d) {
    V = numeric.dot(V, [[1,0,0],[0,1,0],[0,0,-1]])
    console.log("TWIST");
  }
  var U = numeric.dot(V, Wt);

  coords2 = _.map(coords2, function(coord){ return [ coord[0]-Qcx, coord[1]-Qcy, coord[2]-Qcz ] });
  coords2 = numeric.dot(coords2,U);

  coords1 = _.map(coords1, function(coord){ return [ coord[0]-Pcx, coord[1]-Pcy, coord[2]-Pcz ] });

  //Translate coords to stub COM
  //coords1 = _.map(coords1, function(coord){ return [ coord[0]-Pcx, coord[1]-Pcy, coord[2]-Pcz ] });

  //Rotate coords1 onto coords2
  //coords1 = numeric.dot(coords1,rot_matrix);

  //coords2 = _.map(coords2, function(coord){ return [ coord[0]-Qcx, coord[1]-Qcy, coord[2]-Qcz ] });
  //coords1 = _.map(coords1, function(coord){ return [ coord[0]+Qcx, coord[1]+Qcy, coord[2]+Qcz ] });

  //Reapply coordinates, combine PDBS and return
  protein2.change_chain('A', 'B');
  protein1.update_coords_from_matrix(coords1);
  protein2.update_coords_from_matrix(coords2);

  var combined = new Protein();
  combined.chains = protein1.chains.concat(protein2.chains);

  return combined.dump_pdb();
}



///An implementation of the Kabsch algorithm, described in detail
///on wikipedia (http://en.wikipedia.org/wiki/Kabsch_algorithm)
///this function returns U, the optimal rotation matrix, and assumes
///P and Q have already been translated to their center of mass
var kabsch = function(P, Q) {
  var V, S, W, d, U;

  //Find the covariance matrix
  C = numeric.dot(numeric.transpose(P), Q);

  //Calculate the optimal rotation matrix using SVD
  //V,S,W => U,S,V (from numeric.js implementation)
  svd = numeric.svd(C);

  V = svd.U;
  S = svd.S;
  W = numeric.transpose(svd.V);

  //Ignoring d for now assume right handed coordinate system
  d = (numeric.det(V) * numeric.det(W)) < 0.0
  if(d) {
    W = numeric.dot(W, [[1,0,0],[0,1,0],[0,0,-1]])
  }

  //Apply rotation matrix
  U = numeric.dot(V, W);
  return U;
}

var Protein = function() {
  this.chains = [];

  //Get the N,CA,C coords in matrix form for the given residue
  this.get_residue_stub = function(resnum) {
    var residues = [];
    var coords = [];
    _.each(this.chains, function(chain) {
      residues.push($.grep(chain.residues, function(residue) {
        return residue.resnum == resnum;
      })[0]);
    });
    _.each(residues, function(residue) {
      var n = $.grep(residue.atoms, function(atom) {
        return atom.type == "N";
      })[0];
      var ca = $.grep(residue.atoms, function(atom) {
        return atom.type == "CA";
      })[0];
      var c = $.grep(residue.atoms, function(atom) {
        return atom.type == "C";
      })[0];
      coords.push([n.x, n.y, n.z]);
      coords.push([ca.x, ca.y, ca.z]);
      coords.push([c.x, c.y, c.z]);
    })
    return coords;
  }

  this.change_chain = function(old_id, new_id) {
     _.each($.grep(this.chains, function(e){ 
      return e.id === old_id;
    }), function(val) { val.id = new_id });
  }

  this.translate_coords = function(x, y, z) {
    _.each(this.chains, function(chain) { 
      _.each(chain.residues, function(residue) { 
        _.each(residue.atoms, function(atom) {
          atom.x += x;
          atom.y += y;
          atom.z += z;
        })
      })
    });
  }

  this.rotate_coords = function(u) {
    var coords = this.get_coord_matrix()
    var rotated_coords = numeric.dot(coords, u);
    this.update_coords_from_matrix(rotated_coords);
  }

  this.update_coords_from_matrix = function(matrix) {
    var counter = 0;
    _.each(this.chains, function(chain) { 
      _.each(chain.residues, function(residue) { 
        _.each(residue.atoms, function(atom) {
          atom.x = matrix[counter][0];
          atom.y = matrix[counter][1];
          atom.z = matrix[counter][2];
          ++counter;
        })
      })
    });
  }


  ///Get a coordinate matrix from a protein object
  this.get_coord_matrix = function() {
    var coords = [];
    _.each(this.chains, function(chain) { 
      _.each(chain.residues, function(residue) { 
        _.each(residue.atoms, function(atom) {
          coords.push([atom.x, atom.y, atom.z]);
        })
      })
    });
    return coords;
  }

	this.init_from_pdb = function(pdb) {
    var lines = pdb.split('\n');
    for (var i=0; i<lines.length; i+=1) {
      if (lines[i].substr(0,4)=="ATOM" ||
          lines[i].substr(0,6)=="HETATM" ) {
        var x = parseFloat(lines[i].substr(30,7));
        var y = parseFloat(lines[i].substr(38,7));
        var z = parseFloat(lines[i].substr(46,7));
        var chain = lines[i][21].trim();
        var resnum = parseInt(lines[i].substr(22,5).trim());
        var res_type = lines[i].substr(17, 3).trim();
        var atom_type = lines[i].substr(12,4).trim();
        var element = lines[i].substr(76,2).trim();

        var chain_result = $.grep(this.chains, function(e){ return e.id === chain; }); 
        var cur_chain;

        //new chain
        if(chain_result.length === 0) {
          cur_chain = {
            'id': chain,
            'residues': []
          };
          this.chains.push(cur_chain);
        }
        else {
          cur_chain = chain_result[0];
        }

        var residue_result = $.grep(cur_chain.residues, function(e){ return e.resnum === resnum; }); 
        var cur_residue;

        //new chain
        if(residue_result.length === 0) {
          cur_residue = {
            'resnum': resnum,
            'type': res_type,
            'atoms': []
          };
          cur_chain.residues.push(cur_residue);
        }
        else {
          cur_residue = residue_result[0];
        }

        cur_residue.atoms.push({
          'x': x,
          'y': y,
          'z': z,
          'type': atom_type,
          'element': element
        });
      }
    }
  }

  this.dump_pdb = function() {
    var pdb = '';
    var counter = 0;
    _.each(this.chains, function(chain) { 
      _.each(chain.residues, function(residue) { 
        _.each(residue.atoms, function(atom) { 
          ++counter;
          pdb += _.str.sprintf("ATOM  %5d %-4s%1s%-3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f           %-2s%2s\n",
            counter, atom.type, ''/*altLoc*/, residue.type, chain.id, residue.resnum, ''/*icode*/, atom.x, atom.y, atom.z,
            1/*occupancy*/, 1/*temp factor*/, atom.element/*element*/, ''/*charge*/);
        })
      })
    });
    return pdb;
  }
}
var align_and_combine = function(protein1, resnum1, protein2, resnum2) {

  protein2.change_chain('A', 'B');
  var coords1 = protein1.get_coord_matrix();
  var coords2 = protein2.get_coord_matrix();

  //for testing
  protein2.translate_coords(10, 0, 0);

  //Superimpose coordinates 
  var rot_matrix = kabsch(coords1, coords2);
  var rotated_coords = numeric.dot(coords1,rot_matrix);
  protein1.update_coords_from_matrix(rotated_coords);
  //protein2.update_coords_from_matrix(coords2);

  var combined = new Protein();
  combined.chains = protein1.chains.concat(protein2.chains);
  console.log(combined.dump_pdb());

  return combined.dump_pdb();
}



///An implementation of the Kabsch algorithm, described in detail
///on wikipedia (http://en.wikipedia.org/wiki/Kabsch_algorithm)
///this function returns U, the optimal rotation matrix, with P and
///Q having been translated to their centers of mass
var kabsch = function(P, Q) {
  var Pcx, Pcy, Pcz, Qcx, Qcy, Qcz, V, S, W, d, U;

  //Find the center of mass for P and Q, translate coordinates
  Pcx = numeric.sum(_.map(P, function(coord){ return coord[0]; }))/P.length;
  Pcy = numeric.sum(_.map(P, function(coord){ return coord[1]; }))/P.length;
  Pcz = numeric.sum(_.map(P, function(coord){ return coord[2]; }))/P.length;
  P = _.map(P, function(coord){ return [ coord[0]-Pcx, coord[1]-Pcy, coord[2]-Pcz ] });

  Qcx = numeric.sum(_.map(Q, function(coord){ return coord[0]; }))/Q.length;
  Qcy = numeric.sum(_.map(Q, function(coord){ return coord[1]; }))/Q.length;
  Qcz = numeric.sum(_.map(Q, function(coord){ return coord[2]; }))/Q.length;
  Q = _.map(Q, function(coord){ return [ coord[0]-Qcx, coord[1]-Qcy, coord[2]-Qcz ] });

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

  //Apply rotation matrix
  U = numeric.dot(V, W);
  return U;
}

var Protein = function() {
  this.chains = [];

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
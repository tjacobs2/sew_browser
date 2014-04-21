var align_and_combine = function(protein1, resnum1, protein2, resnum2) {
  var coords1 = [];
  var coords2 = [];
  _.each(protein1.chains, function(chain) { 
    _.each(chain.residues, function(residue) { 
      _.each(residue.atoms, function(atom) {
        coords1.push([atom.x, atom.y, atom.z]);
      })
    })
  });

  _.each(protein2.chains, function(chain) { 
    _.each(chain.residues, function(residue) { 
      _.each(residue.atoms, function(atom) {
        coords2.push([atom.x, atom.y, atom.z]);
      })
    })
  });

  protein2.transform_coords(coords2);
  protein2.change_chain('A', 'B');

  //Covariance matrix 

  var P = coords1;
  var Pc = numeric.sum(P)/P.length;
  console.log(Pc);
  var Q = coords2;
  var Qc = numeric.sum(Q)/Q.length;
  console.log(Qc);

  P = numeric.sub(P,Pc);
  Q = numeric.sub(Q,Qc);

  var C = numeric.dot(numeric.transpose(P), Q);
  var svd = numeric.svd(C);

  var V = svd.U;
  var S = svd.S;
  var W = svd.V;
  var d = (numeric.det(V) * numeric.det(W)) < 0.0
  console.log(d);

  //Get rotation matrix
  //V,S,W => U,S,V
  var U = numeric.dot(V, W);
  console.log(U);

  //Rotate P
  var rotated_coords = numeric.dot(P,U);

  console.log(coords1[0]);
  console.log(coords2[0]);
  console.log(rotated_coords[0]);

  var combined = new Protein();
  combined.chains = protein1.chains.concat(protein2.chains);

  return combined.dump_pdb();// + 'TER\n' + protein2.dump_pdb();

}

var Protein = function() {
  this.chains = [];

  this.change_chain = function(old_id, new_id) {
     _.each($.grep(this.chains, function(e){ 
      return e.id === old_id;
    }), function(val) { val.id = new_id });
  }

  this.transform_coords = function(matrix) {
    _.each(this.chains, function(chain) { 
      _.each(chain.residues, function(residue) { 
        _.each(residue.atoms, function(atom) {
          atom.x += 10;
          //atom.y += 1;
          //atom.z += 1;
        })
      })
    });
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
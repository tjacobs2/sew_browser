$(function(){ // on dom ready


//Test out PDB fetching
var glmol = new GLmol('glmol', true);

glmol.defineRepresentation = function() {
   var all = this.getAllAtoms();
   var hetatm = this.removeSolvents(this.getHetatms(all));
   this.colorByAtom(all, {});
   this.colorByChain(all);
   var asu = new THREE.Object3D();
   
   this.drawBondsAsStick(asu, hetatm, this.cylinderRadius, this.cylinderRadius);
   this.drawBondsAsStick(asu, this.getResiduesById(this.getSidechains(this.getChain(all, ['A'])), [58, 87]), this.cylinderRadius, this.cylinderRadius);
   this.drawBondsAsStick(asu, this.getResiduesById(this.getSidechains(this.getChain(all, ['B'])), [63, 92]), this.cylinderRadius, this.cylinderRadius);
   this.drawCartoon(asu, all, this.curveWidth, this.thickness);

   this.drawSymmetryMates2(this.modelGroup, asu, this.protein.biomtMatrices);
   this.modelGroup.add(asu);
};

//Testing out GLMol browser
$.get("http://www.rcsb.org/pdb/files/1yzm.pdb", function(ret1) {
  $.get("http://www.rcsb.org/pdb/files/1yzm.pdb", function(ret2) {
    var protein1 = new Protein();
    var protein2 = new Protein();
    protein1.init_from_pdb(ret1);
    protein2.init_from_pdb(ret2);
    var result = align_and_combine(protein1, 1, protein2, 2);
    $("#glmol_src").val(result);
    glmol.loadMolecule();
  })
});

//retrieve the graph data from the server then render the graph
$.get("/graph", function(ret) {
  var graph_data = ret;
  $('#cy').cytoscape({

    layout: {
      name:'grid'
    },

    style: cytoscape.stylesheet()

      .selector('node')
        .css({
          'content': 'data(id)',
          'text-valign': 'center',
          'color': 'white',
          'text-outline-width': 2,
          'text-outline-color': '#888'
        })
      .selector('edge')
        .css({
          'target-arrow-shape': 'triangle',
          'curve-style': 'haystack'
        })
      .selector(':selected')
        .css({
          'background-color': 'black',
          'line-color': 'black',
          'target-arrow-color': 'black',
          'source-arrow-color': 'black'
        })
      .selector('.faded')
        .css({
          'opacity': 0.25,
          'text-opacity': 0
        }),
    
    ready: function(){
      window.cy = this;
      
      // giddy up...
      cy.load(graph_data);
      cy.boxSelectionEnabled(false);
      cy.userPanningEnabled(true);
      cy.userZoomingEnabled(true);
      
      cy.elements().unselectify();

      cy.on('tap', 'edge', function(e){
        var edge = e.cyTarget; 
        //var node_1 = edge.source();
        var neighborhood = edge.target().closedNeighborhood().add(edge.source().closedNeighborhood());

        //console.log(node_1);
        // var neighborhood = node.neighborhood().add(node);
        
        //cy.elements().addClass('faded');
        neighborhood.removeClass('faded');
      });
      
      //When you tap an edge, show its neighborhood and load the corresponding pdb
      cy.on('tap', 'node', function(e){
        var node = e.cyTarget; 
        var neighborhood = node.neighborhood().add(node);
        
        cy.elements().addClass('faded');
        neighborhood.removeClass('faded');

        var pdb_path = node.data().pdb_code.split('/');
        var pdb_code = pdb_path[pdb_path.length - 1].split('_')[0];
        console.log("Fetching "+pdb_code+" from PDB");

        $.get("http://www.rcsb.org/pdb/files/"+pdb_code+".pdb", function(ret) {
          var protein = new Protein();
          protein.init_from_pdb(ret);
          var pdb = protein.dump_pdb();
          $("#glmol_src").val(pdb);
          //console.log(pdb);
          glmol.loadMolecule();
        });

        //center a node on click
        //cy.center(node);
      });
      
      cy.on('tap', function(e){
        if( e.cyTarget === cy ){
          cy.elements().removeClass('faded');
        }
      });
    }
  });
});

}); // on dom ready
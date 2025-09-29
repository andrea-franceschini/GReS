classdef MeshGlueJumpStabilization < MeshGlue
  % Mesh glue class implementing pressure jump stabilization

  properties (Access = private)
    scale
  end

  methods (Access = public)
    function obj = MeshGlueJumpStabilization(id,inputStruct,domains)
      obj@MeshGlue(id,inputStruct,domains);
      if isfield(inputStruct.Stabilization,"scaleAttribute")
        obj.scale = inputStruct.Stabilization.scaleAttribute;
      else 
        obj.scale = 1.0;
      end
      obj.multiplierType = 'P0';
    end

    function computeMat(obj,dt)
      computeMat@MeshGlue(obj,dt);

      % map local mortar matrices to global indices
      obj.Jmult = -computeStabilizationMatrix(obj,obj.physic);
    end

end



  methods(Access = private)
    function out = isStabMatrixReady(obj)
      out = all(...
        [~cellfun(@isempty, obj.Jslave), ~cellfun(@isempty, obj.Jmaster)]);
    end

    function stabMat = computeStabilizationMatrix(obj,fld)

      % get number of components of input field
      nc = obj.dofm(1).getDoFperEnt(fld);

      % initialize matrix estimating number of entries
      % number of internal slave elements
      nes = sum(all(obj.mesh.e2f{2},2));
      nEntries = 2*nc*nes; % each cell should contribute at least two times
      [id1,id2,vals] = deal(zeros(nEntries,1));

      c = 0;

      % get list of internal master edges
      inEdgeMaster = find(all(obj.mesh.e2f{1},2));

      for ieM = inEdgeMaster'
        % get master faces sharing internal edge ie
        fM = obj.mesh.e2f{1}(ieM,:);
        assert(numel(fM)==2,['Unexpected number of connected faces for' ...
          'master edge %i. Expected 2.'], ieM);

        % get slave faces sharing support with master faces
        fS = unique([find(obj.mesh.elemConnectivity(fM(1),:)),...
          find(obj.mesh.elemConnectivity(fM(2),:))]);

        if numel(fS) < 2
          continue
        end

        % get internal edges of slave faces

        eS = unique(obj.mesh.f2e{2}(fS,:));
        id = all(ismember(obj.mesh.e2f{2}(eS,:),fS),2);
        ieS = eS(id);

        % get active macroelement nodes
        nM = obj.mesh.e2n{1}(ieM,:);
        nS = unique(obj.mesh.e2n{2}(eS,:));

        % compute local schur complement approximation
        Stmp = computeSchurLocal(obj,nM,nS,fS,fld);
        S = zeros(nc,1);
        for i = 1:nc
          S(i) = mean(Stmp(i:nc:end));
        end

        % apply scaling due to relative grid size
        Am = mean(obj.mesh.msh(1).surfaceArea(fM));
        As = mean(obj.mesh.msh(2).surfaceArea(fS));
        S = (Am/As)*S;

        % assemble stabilization matrix component
        for iesLoc = ieS'
          f = obj.mesh.e2f{2}(iesLoc,:);
          id1(c+1:c+nc) = dofId(f(1),nc);
          id2(c+1:c+nc) = dofId(f(2),nc);
          vals(c+1:c+nc) = S;
          c = c+nc;
        end
      end

      id1 = id1(1:c); id2 = id2(1:c); vals = vals(1:c);
      % assemble sparse matrix
      nmult = nc*obj.mesh.msh(2).nSurfaces;
      stabMat = sparse(id1,id1,vals,nmult,nmult)+...
        sparse(id1,id2,-vals,nmult,nmult)+...
        sparse(id2,id2,vals,nmult,nmult);
      stabMat = stabMat + stabMat' - diag(diag(stabMat));
    end

    function S = computeSchurLocal(obj,nm,ns,fs,field)
      % compute approximate schur complement for local nonconforming
      % patch of element
      % input: nm/ns local master/slave node indices
      % fs: local slave faces indices

      % get slave and master dof to access jacobian
      fldM = getFieldId(obj.dofm(1),field);
      fldS = getFieldId(obj.dofm(2),field);
      dofS = obj.dofm(2).getLocalDoF(obj.mesh.local2glob{2}(ns),fldS);
      dofM = obj.dofm(1).getLocalDoF(obj.mesh.local2glob{1}(nm),fldM);
      nc = obj.dofm(1).getDoFperEnt(field);


      % get local mortar matrices
      Dloc = obj.Jslave(dofId(fs,nc),dofS);
      Mloc = obj.Jmaster(dofId(fs,nc),dofM);
      V = [Dloc, Mloc];              % minus sign!
      %V = Discretizer.expandMat(V,nc);

      % get local jacobian
      Km = getSolver(obj.solvers(1),field).J(dofM,dofM);
      Ks = getSolver(obj.solvers(2),field).J(dofS,dofS);
      Kloc = diag([1./diag(Ks);1./diag(Km)]);

      S = obj.scale*V*(Kloc*V');  % compute Schur complement
    end

  end
end


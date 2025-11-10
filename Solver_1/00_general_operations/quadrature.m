function [quadPoints,quadWeights] = quadrature(elemType,numberOfPoints)
% Calculates the position and weight of the gauss points in the reference element
% INPUT
%   elemType         0=line, 1=square, 2=triangle
%   numberOfPoints   number of quadrature points
% OUTPUT
%   position and weigth of gauss points in the reference element
%

switch(elemType)
   case 0 % 1D lines
      switch(numberOfPoints)
         case 2 % degree 2
            pos1 = 1/sqrt(3);
            quadPoints = [-pos1 pos1 ];
            quadWeights = [1 1];
          case 3 % degree 2
              pos1 = sqrt(3/5);
              quadPoints = [-pos1 0 pos1];
              quadWeights = [5/9 8/9 5/9];
          case 4
              a1 = (3/7)+(2/7)*sqrt(6/5);
              b1 = (3/7)-(2/7)*sqrt(6/5);
              pos1 = sqrt(a1);
              pos2 = sqrt(b1);
              quadPoints = [-pos1 -pos2 pos2 pos1];
              pg1 = (18-sqrt(30))/36;
              pg2 = (18+sqrt(30))/36;
              pg3 = pg2;
              pg4 = pg1;
              quadWeights = [pg1 pg2 pg3 pg4];
          case 5
              a1 = 5+2*sqrt(10/7); 
              b1 = 5-2*sqrt(10/7);
              c1 = 0;
              pos1 = (1/3)*sqrt(a1);
              pos2 = (1/3)*sqrt(b1);
              pos3 = c1;
              quadPoints = [-pos1 -pos2 pos3 pos2 pos1];
              pg1 = (322-13*sqrt(70))/900;
              pg2 = (322+13*sqrt(70))/900;
              pg3 = 128/225;
              quadWeights = [pg1 pg2 pg3 pg2 pg1];              
          case 9
              pos1 =  0.96816023950762608983557620290367;
              pos2 =  0.83603110732663579429942978806973;
              pos3 =  0.61337143270059039730870203934147;
              pos4 =  0.32425342340380892903853801464334;
              pos5 = 0;
              quadPoints = [-pos1 -pos2 -pos3 -pos4 pos5 pos4 pos3 pos2 pos1];
              pg1 = 0.081274388361574411971892158110524;
              pg2 = 0.18064816069485740405847203124291;
              pg3 = 0.26061069640293546231874286941863;
              pg4 = 0.31234707704000284006863040658444;
              pg5 = 32768/99225;
              quadWeights = [pg1 pg2 pg3 pg4 pg5 pg4 pg3 pg2 pg1];
         otherwise
            error('Quadrature for linear elements with %i points not implemented', numberOfPoints)
      end
      
   case 1 % squares
      switch(numberOfPoints)
         case 4 % degree 1
            pos1 = 1/sqrt(3);
            quadPoints = [-pos1 -pos1;  pos1 -pos1; 
                           pos1  pos1; -pos1  pos1];
            quadWeights = [1 1 1 1];
         case 9  % degree 2
            pos1 = sqrt(3/5);
            quadPoints = [-pos1 -pos1; 0 -pos1; pos1 -pos1; 
                          -pos1     0; 0     0; pos1     0; 
                          -pos1  pos1; 0  pos1; pos1  pos1];
            pg1 = 5/9;
            pg2 = 8/9;
            pg3 = pg1;
            quadWeights= [pg1*pg1, pg2*pg1, pg3*pg1, pg1*pg2, pg2*pg2, ...
               pg3*pg2, pg1*pg3, pg2*pg3, pg3*pg3];
        
          case 16  % degree 3
              a1 = (3/7)+(2/7)*sqrt(6/5);
              b1 = (3/7)-(2/7)*sqrt(6/5);
              pos1 = sqrt(a1);
              pos2 = sqrt(b1);
              quadPoints = [-pos1 -pos1; -pos2 -pos1; pos2 -pos1; pos1 -pos1;
                  -pos1 -pos2; -pos2 -pos2; pos2 -pos2; pos1 -pos2; 
                  -pos1 pos2; -pos2 pos2; pos2 pos2; pos1 pos2; 
                  -pos1 pos1; -pos2 pos1; pos2 pos1; pos1 pos1];

              pg1 = (18-sqrt(30))/36;
              pg2 = (18+sqrt(30))/36;
              pg3 = pg2;
              pg4 = pg1;
              quadWeights= [pg1*pg1, pg2*pg1, pg3*pg1, pg4*pg1, ...
                  pg1*pg2, pg2*pg2, pg3*pg2, pg4*pg2, ...
                  pg1*pg3, pg2*pg3, pg3*pg3, pg4*pg3, ...
                  pg1*pg4, pg2*pg4, pg3*pg4, pg4*pg4];
         case 25  % degree 5
             a1 = 5+2*sqrt(10/7); 
             b1 = 5-2*sqrt(10/7);
             c1 = 0;
              pos1 = (1/3)*sqrt(a1);
              pos2 = (1/3)*sqrt(b1);
              pos3 = c1;
              quadPoints = [-pos1 -pos1; -pos2 -pos1; pos3 -pos1; pos2 -pos1; pos1 -pos1;
                  -pos1 -pos2; -pos2 -pos2; pos3 -pos2; pos2 -pos2; pos1 -pos2;
                  -pos1 pos3; -pos2 pos3; pos3 pos3; pos2 pos3; pos1 pos3;
                  -pos1 pos2; -pos2 pos2; pos3 pos2; pos2 pos2; pos1 pos2;
                  -pos1 pos1; -pos2 pos1; pos3 pos1; pos2 pos1; pos1 pos1];

              pg1 = (322-13*sqrt(70))/900;
              pg2 = (322+13*sqrt(70))/900;
              pg3 = 128/225;
              pg4 = pg2;
              pg5 = pg1;
              quadWeights= [pg1*pg1,pg2*pg1,pg3*pg1,pg4*pg1,pg5*pg1,...
                  pg1*pg2,pg2*pg2,pg3*pg2,pg4*pg2,pg5*pg2,...
                  pg1*pg3,pg2*pg3,pg3*pg3,pg4*pg3,pg5*pg3,...
                  pg1*pg4,pg2*pg4,pg3*pg4,pg4*pg4,pg5*pg4,...
                  pg1*pg5,pg2*pg5,pg3*pg5,pg4*pg5,pg5*pg5];
          case 81 % degree 9
              pos1 =  0.96816023950762608983557620290367;
              pos2 =  0.83603110732663579429942978806973;
              pos3 =  0.61337143270059039730870203934147;
              pos4 =  0.32425342340380892903853801464334;
              pos5 = 0;
              quadPoints = [-pos1 -pos1; -pos2 -pos1; -pos3 -pos1; -pos4 -pos1; pos5 -pos1; pos4 -pos1; pos3 -pos1; pos2 -pos1; pos1 -pos1;
                  -pos1 -pos2; -pos2 -pos2; -pos3 -pos2; -pos4 -pos2; pos5 -pos2; pos4 -pos2; pos3 -pos2; pos2 -pos2; pos1 -pos2;
                  -pos1 -pos3; -pos2 -pos3; -pos3 -pos3; -pos4 -pos3; pos5 -pos3; pos4 -pos3; pos3 -pos3; pos2 -pos3; pos1 -pos3;
                  -pos1 -pos4; -pos2 -pos4; -pos3 -pos4; -pos4 -pos4; pos5 -pos4; pos4 -pos4; pos3 -pos4; pos2 -pos4; pos1 -pos4;
                  -pos1 pos5; -pos2 pos5; -pos3 pos5; -pos4 pos5; pos5 pos5; pos4 pos5; pos3 pos5; pos2 pos5; pos1 pos5;
                  -pos1 pos4; -pos2 pos4; -pos3 pos4; -pos4 pos4; pos5 pos4; pos4 pos4; pos3 pos4; pos2 pos4; pos1 pos4;
                  -pos1 pos3; -pos2 pos3; -pos3 pos3; -pos4 pos3; pos5 pos3; pos4 pos3; pos3 pos3; pos2 pos3; pos1 pos3;
                  -pos1 pos2; -pos2 pos2; -pos3 pos2; -pos4 pos2; pos5 pos2; pos4 pos2; pos3 pos2; pos2 pos2; pos1 pos2;
                  -pos1 pos1; -pos2 pos1; -pos3 pos1; -pos4 pos1; pos5 pos1; pos4 pos1; pos3 pos1; pos2 pos1; pos1 pos1];
              pg1 = 0.081274388361574411971892158110524;
              pg2 = 0.18064816069485740405847203124291;
              pg3 = 0.26061069640293546231874286941863;
              pg4 = 0.31234707704000284006863040658444;
              pg5 = 32768/99225;
              pg6 = pg4;
              pg7 = pg3;
              pg8 = pg2;
              pg9 = pg1;
              quadWeights = [pg1*pg1, pg2*pg1, pg3*pg1, pg4*pg1, pg5*pg1, pg6*pg1, pg7*pg1, pg8*pg1, pg9*pg1,...
                                        pg1*pg2, pg2*pg2, pg3*pg2, pg4*pg2, pg5*pg2, pg6*pg2, pg7*pg2, pg8*pg2, pg9*pg2,...
                                      pg1*pg3, pg2*pg3, pg3*pg3, pg4*pg3, pg5*pg3, pg6*pg3, pg7*pg3, pg8*pg3, pg9*pg3,...
                                      pg1*pg4, pg2*pg4, pg3*pg4, pg4*pg4, pg5*pg4, pg6*pg4, pg7*pg4, pg8*pg4, pg9*pg4,...
                                      pg1*pg5, pg2*pg5, pg3*pg5, pg4*pg5, pg5*pg5, pg6*pg5, pg7*pg5, pg8*pg5, pg9*pg5,...
                                      pg1*pg6, pg2*pg6, pg3*pg6, pg4*pg6, pg5*pg6, pg6*pg6, pg7*pg6, pg8*pg6, pg9*pg6,...
                                      pg1*pg7, pg2*pg7, pg3*pg7, pg4*pg7, pg5*pg7, pg6*pg7, pg7*pg7, pg8*pg7, pg9*pg7,...
                                      pg1*pg8, pg2*pg8, pg3*pg8, pg4*pg8, pg5*pg8, pg6*pg8, pg7*pg8, pg8*pg8, pg9*pg8,...
                                      pg1*pg9, pg2*pg9, pg3*pg9, pg4*pg9, pg5*pg9, pg6*pg9, pg7*pg9, pg8*pg9, pg9*pg9];

         otherwise
            error('Quadrature for squares with %i points not implemented', numberOfPoints)
      end
      
   case 2 % triangles
      switch(numberOfPoints)
         case 1  % degree 1
            pos1 = 1/3;
            quadPoints = [pos1  pos1];
            quadWeights = 1/2;
         case 3  % degree 2
            pos1 = 1/2;
            quadPoints = [pos1 pos1; 0 pos1; pos1 0];
            pes1 = 1/6;
            quadWeights = [pes1 pes1 pes1];
         case 4  % degree 3
            quadPoints = [1/3 1/3; 0.6 0.2; 0.2 0.6; 0.2 0.2];
            quadWeights = [-27/96   25/96   25/96   25/96];
         case 7   % degree 5
            a = 0.101286507323456338800987361915123;
            b = 0.470142064105115089770441209513447;
            P1 = 0.0629695902724135762978419727500906;
            P2 = 0.0661970763942530903688246939165759;
            quadPoints = [a a; a 1-2*a; 1-2*a a; b b; b 1-2*b; 1-2*b b; ...
               1/3 1/3];
            quadWeights = [P1, P1, P1, P2, P2, P2, 0.1125];
         case 12   % degree 6
            a = 0.0630890144915022283403316028708191;
            b = 0.249286745170910421291638553107019;
            c = 0.0531450498448169473532496716313981;
            d = 0.310352451033784405416607733956552;
            P1 = 0.0254224531851034084604684045534344;
            P2 = 0.0583931378631896830126448056927897;
            P3 = 0.0414255378091867875967767282102212;
            quadPoints = [a a; a 1-2*a; 1-2*a a; b b; b 1-2*b; 1-2*b b; ...
               c d; c 1-c-d; 1-c-d c; 1-c-d d; d 1-c-d; d c];
            quadWeights = [P1, P1, P1, P2, P2, P2, P3, P3, P3, P3, P3, P3];
         case 16  % degree 8
            a = 0.170569307751760206622293501491464;
            b = 0.0505472283170309754584235505965989;
            c = 0.459292588292723156028815514494169;
            d = 0.728492392955404281241000379176061;
            e = 0.263112829634638113421785786284643;
            f = 1/3;
            P1 = 0.0516086852673591251408957751460645;
            P2 = 0.0162292488115990401554629641708902;
            P3 = 0.0475458171336423123969480521942921;
            P4 = 0.0136151570872174971324223450369544;
            P5 = 0.0721578038388935841255455552445323;
            quadPoints = [a a; a 1-2*a; 1-2*a a; b b; b 1-2*b; 1-2*b b;...
               c c; c 1-2*c; 1-2*c c; d e; d 1-d-e; 1-d-e d; 1-d-e e; ...
               e 1-d-e; e d; f f];
            quadWeights = [P1*[1 1 1], P2*[1 1 1], P3*[1 1 1], ...
               P4*ones(1,6), P5];
         otherwise
            error('Quadrature for triangles with %i points not implemented', numberOfPoints)
      end
   otherwise
      error('Unknown element type: %i \nValid options are 0=1D lines, 1=squares, 2=triangles', elemType)
end

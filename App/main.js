var Objetos = [];
var camera;
var iluminacao;
let I;
let telaX = window.screen.availWidth/1.5;
let telaY = window.screen.availHeight;
var ZbufferX = [];
var Zbuffer = [];
//var bufferCorX = [];
var bufferCor = [];
let Cd;
let Cs;
let Ca;
let CIlum;
let qtdpint = 0;

function Inicializar(){
    InicializadorZbuffer();     
    prepararCamera();
    mudarCoordenadasVista();
    normais_Triangulos_Vertices();
    normalizar_Normais();
    projecaoPespectiva();
    mudarCoordenadasTela();
    ordenarVertices();
    desenharTriangulos();
    draw();

    console.log(Zbuffer);
    console.log(bufferCor);
    console.log(Objetos);
    console.log(qtdpint);
    
}

function InicializadorZbuffer(){
    for (let i = 0; i < telaY; i++) {
        for (let j = 0; j < telaX; j++) {
            ZbufferX[j] = Infinity;
            //bufferCorX[j] = new Vetor(0,0,0);
        }
        Zbuffer[i] = ZbufferX;
        //bufferCor[i] = bufferCorX;
    }
}

//Função para criar o vetor U e a matriz I
function prepararCamera(){

    camera.V_orto = OrtogonalizarVetor(camera.V, camera.N);
    camera.V_norm = normalizarVetor(camera.V_orto);
    camera.N_norm = normalizarVetor(camera.N);
    camera.U = produtoVetorial(camera.V_norm, camera.N_norm);

    I=matrizVetores(camera.U,camera.V_norm,camera.N_norm);

}

//Função para mudar as coordenadas de mundo para de vista
function mudarCoordenadasVista(){
    //Mudar coordenadas dos pontos
    for (let i = 0; i < Objetos.length; i++) {
        for (let j = 0; j < Objetos[i].qtd_pontos; j++) {
            let aux = Objetos[i].pontos[j];
            let aux_vet = SubtracaoPontos(camera.C,aux);
            Objetos[i].pontos[j] = MultMatrizVetor(I,aux_vet);
        }
    }
    //Mudar coordenadas da fonte de luz
    let auxIlu_vet = SubtracaoPontos(camera.C,iluminacao.Pl);
    iluminacao.Pl = MultMatrizVetor(I,auxIlu_vet);

}

//Função para calcular as normais dos triangulos e das suas vertices
function normais_Triangulos_Vertices(){
    Objetos[Objetos.length-1].NormalTriangulo();
    //Objetos[Objetos.length-1].NormalVertices();
}

//Função para normalizar as normais dos triangulos e vertices
function normalizar_Normais(){
    //Normaliza as normas dos triangulos
    for (let i = 0; i < Objetos.length; i++) {
        for (let j = 0; j < Objetos[i].triangulos.length; j++) {
            let aux = Objetos[i].triangulos[j].normal;
            Objetos[i].triangulos[j].normal = normalizarVetor(aux);
        }
    }
    //Normaliza as normas dos vertices
    for (let i = 0; i < Objetos.length; i++) {
        for (let j = 0; j < Objetos[i].pontos.length; j++) {
            let aux = Objetos[i].pontos[j].normal;
            Objetos[i].pontos[j].normal = normalizarVetor(aux);
        }
    }
    
}

//Função para converter os vertices de coordenadas de vista para perspectiva
function projecaoPespectiva(){
    for (let i = 0; i < Objetos.length; i++) {
        for (let j = 0; j < Objetos[i].pontos.length; j++) {
            Objetos[i].pontos_tela[j] = new Ponto();
            Objetos[i].pontos_tela[j].x = ((Objetos[i].pontos[j].x/Objetos[i].pontos[j].z)*(camera.d/camera.hx));
            Objetos[i].pontos_tela[j].y = ((Objetos[i].pontos[j].y/Objetos[i].pontos[j].z)*(camera.d/camera.hy));
            Objetos[i].pontos_tela[j].z = 0;
        }
    }
}

//Função para converter para as coordenadas de tela
function mudarCoordenadasTela(){
        
    for (let i = 0; i < Objetos.length; i++) {
        for (let j = 0; j < Objetos[i].pontos.length; j++) {
            
            let ptX = Objetos[i].pontos_tela[j].x;
            let ptY = Objetos[i].pontos_tela[j].y;
           
            Objetos[i].pontos_tela[j].x = parseInt(((ptX+1)*telaX)/2);
            Objetos[i].pontos_tela[j].y = parseInt(((1-ptY)*telaY)/2);

        }
    }
}

//Função para Ordenar os vertices dos triangulos pelo menor Y
function ordenarVertices(){
    for (let i = 0; i < Objetos.length; i++) {
        for (let j = 0; j < Objetos[i].triangulos.length; j++) {

            let aux1=Objetos[i].triangulos[j].vertice1;
            let aux2=Objetos[i].triangulos[j].vertice2;
            let aux3=Objetos[i].triangulos[j].vertice3;

            if(Objetos[i].pontos[aux1].y<Objetos[i].pontos[aux2].y){
                Objetos[i].triangulos[j].vertice1 = aux2;
                Objetos[i].triangulos[j].vertice2 = aux1;
                aux1=Objetos[i].triangulos[j].vertice1;
                aux2=Objetos[i].triangulos[j].vertice2;

            }else if(Objetos[i].pontos[aux1].y==Objetos[i].pontos[aux2].y){
                if(Objetos[i].pontos[aux1].x<Objetos[i].pontos[aux2].x){
                    Objetos[i].triangulos[j].vertice1 = aux2;
                    Objetos[i].triangulos[j].vertice2 = aux1;
                    aux1=Objetos[i].triangulos[j].vertice1;
                    aux2=Objetos[i].triangulos[j].vertice2;
                }
            }

            if(Objetos[i].pontos[aux1].y<Objetos[i].pontos[aux3].y){
                Objetos[i].triangulos[j].vertice1 = aux3;
                Objetos[i].triangulos[j].vertice3 = aux1;
                aux1=Objetos[i].triangulos[j].vertice1;
                aux3=Objetos[i].triangulos[j].vertice3;

            }else if(Objetos[i].pontos[aux1].y==Objetos[i].pontos[aux2].y){
                if(Objetos[i].pontos[aux1].x<Objetos[i].pontos[aux2].x){
                    Objetos[i].triangulos[j].vertice1 = aux3;
                    Objetos[i].triangulos[j].vertice3 = aux1;
                    aux1=Objetos[i].triangulos[j].vertice1;
                    aux3=Objetos[i].triangulos[j].vertice3;
                }
            }

            if(Objetos[i].pontos[aux2].y<Objetos[i].pontos[aux3].y){
                Objetos[i].triangulos[j].vertice2 = aux3;
                Objetos[i].triangulos[j].vertice3 = aux2;
                aux2=Objetos[i].triangulos[j].vertice2;
                aux3=Objetos[i].triangulos[j].vertice3;

            }else if(Objetos[i].pontos[aux2].y==Objetos[i].pontos[aux3].y){
                if(Objetos[i].pontos[aux2].x<Objetos[i].pontos[aux3].x){
                    Objetos[i].triangulos[j].vertice2 = aux3;
                    Objetos[i].triangulos[j].vertice3 = aux2;
                    aux2=Objetos[i].triangulos[j].vertice2;
                    aux3=Objetos[i].triangulos[j].vertice3;
                }
            }
        }
    }
}

function desenharTrianguloYIgualBaixo(objetoIndice,vt1,vt2,vt3,bo){
    let v1 = Objetos[objetoIndice].pontos_tela[vt1];
    let v2 = Objetos[objetoIndice].pontos_tela[vt2];
    let v3 = vt3;

    if(bo){
        v3 = Objetos[objetoIndice].pontos_tela[vt3];
    }

    let Xmin = v1.x;
    let Xmax = v1.x;

    let Mmin = parseInt(parseFloat(v2.y-v1.y)/parseFloat(v2.x-v1.x));
    let Mmax = parseInt(parseFloat(v3.y-v1.y)/parseFloat(v3.x-v1.x));

    for (let Yscan = v1.y; Yscan < v3.y; Yscan++) {
        if(bo){
            //desenhador(Xmin,Xmax,Yscan,vt1,vt2,vt3,objetoIndice)
        }else{
            desenhadorBaixo(Xmin,Xmax,Yscan,vt1,vt2,vt3,objetoIndice);
        }
        Xmin += parseInt(1/parseFloat(Mmin));
        Xmax += parseInt(1/parseFloat(Mmax));

    }

}

function desenharTrianguloYIgualCima(objetoIndice,vt1,vt2,vt3,bo){
    let v1 = Objetos[objetoIndice].pontos_tela[vt1];
    let v2 = vt2;
    let v3 = Objetos[objetoIndice].pontos_tela[vt3];

    if(bo){
        v2 = Objetos[objetoIndice].pontos_tela[vt2];
    }

    let Xmin = v3.x;
    let Xmax = v3.x;

    let Mmin = parseInt(parseFloat(v3.y-v1.y)/parseFloat(v3.x-v1.x));
    let Mmax = parseInt(parseFloat(v3.y-v2.y)/parseFloat(v3.x-v2.x));

    for (let Yscan = v3.y; Yscan > v1.y; Yscan--) {
        if(bo){
            //desenhador(Xmin,Xmax,Yscan,vt1,vt2,vt3,objetoIndice);
        }else{
            desenhadorCima(Xmin,Xmax,Yscan,vt1,vt2,vt3,objetoIndice);
        }
        Xmin -= parseInt(1/parseFloat(Mmin));
        Xmax -= parseInt(1/parseFloat(Mmax));

    }
}

function desenharTriangulos(){
    for (let i = 0; i < Objetos.length; i++) {
        for (let j = 0; j < Objetos[i].triangulos.length; j++) {

            let v1 = Objetos[i].triangulos[j].vertice1;
            let v2 = Objetos[i].triangulos[j].vertice2;
            let v3 = Objetos[i].triangulos[j].vertice3;

            let p1 = Objetos[i].pontos[v1];
            let p2 = Objetos[i].pontos[v2];
            let p3 = Objetos[i].pontos[v3];

            let Xmim = Objetos[i].pontos_tela[v1].x;
            let Xmax = Objetos[i].pontos_tela[v1].x;
            let Ymim = Objetos[i].pontos_tela[v1].y;
            let Ymax = Objetos[i].pontos_tela[v3].y;

            let Mv12 = (parseFloat(p2.y-p1.y)/parseFloat(p2.x-p1.x));
            let Mv13 = (parseFloat(p3.y-p1.y)/parseFloat(p3.x-p1.x));
            let Mv23 = (parseFloat(p3.y-p2.y)/parseFloat(p3.x-p2.x));
            
            let verificado=true;

            if(p1.y == p2.y){
                Xmax = p2.x;
                Mv12 = Mv23;

                verificado = false;
            }else if(p1.y == p3.y){
                Xmax = p3.x;
                Mv13 = Mv23;
                verificado = false;
            }

            scanline(Xmim, Xmax, Ymim, Ymax, v1, v2, v3, Mv12, Mv13, Mv23, i, verificado);

            /*let tri = Objetos[i].triangulos[j];

            if(Objetos[i].pontos_tela[tri.vertice1].y==Objetos[i].pontos_tela[tri.vertice2].y){
                desenharTrianguloYIgualCima(i,tri.vertice1,tri.vertice2,tri.vertice3,true);

            }else if(Objetos[i].pontos_tela[tri.vertice2].y == Objetos[i].pontos_tela[tri.vertice3].y){
                desenharTrianguloYIgualBaixo(i,tri.vertice1,tri.vertice2,tri.vertice3,true);

            }else{
                let vertice_aux = new Ponto();
                vertice_aux.x = parseInt(Objetos[i].pontos_tela[tri.vertice1].x + ((parseFloat(Objetos[i].pontos_tela[tri.vertice2].y - Objetos[i].pontos_tela[tri.vertice1].y)/parseFloat(Objetos[i].pontos_tela[tri.vertice3].y-Objetos[i].pontos_tela[tri.vertice1].y))*(Objetos[i].pontos_tela[tri.vertice3].x-Objetos[i].pontos_tela[tri.vertice1].x)));
                vertice_aux.y = Objetos[i].pontos_tela[tri.vertice2].y;
                vertice_aux.z = ( Math.abs((parseFloat(Objetos[i].pontos[tri.vertice1].z - Objetos[i].pontos[tri.vertice3].z))/2) + Objetos[i].pontos[tri.vertice1].z);
                vertice_aux.normal = Objetos[i].triangulos[j].normal;
                desenharTrianguloYIgualBaixo(i,tri.vertice1,tri.vertice2, vertice_aux,false);
                desenharTrianguloYIgualCima(i,tri.vertice2,vertice_aux,tri.vertice3,false);
            }*/
            
        }
    }
}

function scanline(Xmin,Xmax,Ymim,Ymax,v1,v2,v3,Mv12, Mv13, Mv23, objetoIndice, verificado){
    
    let pt1 = Objetos[objetoIndice].pontos_tela[v1];
    let pt2 = Objetos[objetoIndice].pontos_tela[v2];
    let pt3 = Objetos[objetoIndice].pontos_tela[v3];

    let p1 = Objetos[objetoIndice].pontos[v1];
    let p2 = Objetos[objetoIndice].pontos[v2];
    let p3 = Objetos[objetoIndice].pontos[v3];

    let bari = [];
    let p_linha = new Ponto();
    
    for (let y = Ymim; y <= Ymax; y++) {
        
        
        for (let x = Xmin; x <= Xmax; x++) {
            x = Math.floor(x);
            y = Math.floor(y);
            
            if(x>0 && x<telaX && y>0 && y<telaY){
                let p = new Ponto(0,0,0);
                p.x = x;
                p.y = y;

                bari = acharBaricentro(p,pt1,pt2,pt3);

                p_linha = SomaPontos(SomaPontos(multiplicarPontoPorEscalar(bari[0],p1),multiplicarPontoPorEscalar(bari[1],p2)),multiplicarPontoPorEscalar(bari[2],p3));
                //console.log(p);
                //console.log(p_linha);
                
                //console.log(Zbuffer);
                
                if( Zbuffer.length > p.x && Zbuffer[p.x].length > p.y && p_linha.z <= Zbuffer[p.x][p.y]){
                    qtdpint++;
                    Zbuffer[p.x][p.y] = p_linha.z;
                    let zero = new Ponto(0,0,0);
                    let N = SomaVetores(SomaVetores(multiplicarVetorPorEscalar(bari[0],p1.normal),multiplicarVetorPorEscalar(bari[1],p2.normal)),multiplicarVetorPorEscalar(bari[2],p3.normal));
                    let V = SubtracaoPontos(zero, p_linha);
                    let L = SubtracaoVetorPonto(iluminacao.Pl, p_linha);
                    let R = SubtracaoVetores( multiplicarVetorPorEscalar((2*produtoInterno(L,N)),N) ,L);/* */

                    N = normalizarVetor(N);
                    V = normalizarVetor(V);
                    L = normalizarVetor(L);
                    R = normalizarVetor(R);/* */

                    if(produtoInterno(N,V)<0){
                        N = multiplicarVetorPorEscalar(-1,N);
                    }
        
                    Ca = multiplicarVetorPorEscalar(iluminacao.ka, iluminacao.Ia);
                    Cd = MultiplicacaoComponenteComponente(multiplicarVetorPorEscalar((iluminacao.kd * produtoInterno(N,L)),iluminacao.Od),iluminacao.Il);
                    Cs = multiplicarVetorPorEscalar((iluminacao.ks*(Math.pow(produtoInterno(R,V),iluminacao.n))),iluminacao.Il);

                    if(produtoInterno(N,L)<0){
                        Cd = new Vetor(0,0,0);
                        Cs = new Vetor(0,0,0);
                    }else{
                        if(produtoInterno(R,V)<0){
                            Cs = new Vetor(0,0,0);
                        }
                    }
                    
                    CIlum = SomaVetores(SomaVetores(Ca,Cd),Cs);

                    if(CIlum.x > 255){
                        CIlum.x = 255;
                    }
                    if(CIlum.y > 255){
                        CIlum.y = 255;
                    }
                    if(CIlum.z > 255){
                        CIlum.z = 255;
                    }
                    
                    let cor = {
                        r: CIlum.x,
                        g: CIlum.y,
                        b: CIlum.z
                    }

                    //console.log(CIlum);
                    
                    //bufferCor[p.x][p.y] = CIlum;
                    bufferCor.push({x:p.x,y:p.y, cor:cor});
                    //draw();
                    //console.log(Ca,Cd,Cs);
                }

            }


        }

        if(verificado && (y==pt2.y)||(y==pt3.y)){
            if (y==pt2.y) {
                Mv12 = Mv23;
            }else{
                Mv13 = Mv23;
            }
            verificado = false;
        }

        if (Mv12 != Infinity && Mv12 != -Infinity) {
            Xmin += 1/Mv12;
        }

        if (Mv13 != Infinity && Mv13 != -Infinity) {
            Xmax += 1/Mv13;
        }

    }


}

function desenhador(Xmin,Xmax,Yscan,vt1,vt2,vt3,objetoIndice){
    let v1 = Objetos[objetoIndice].pontos_tela[vt1];
    let v2 = Objetos[objetoIndice].pontos_tela[vt2];
    let v3 = Objetos[objetoIndice].pontos_tela[vt3];

    let vv1 = Objetos[objetoIndice].pontos[vt1];
    let vv2 = Objetos[objetoIndice].pontos[vt2];
    let vv3 = Objetos[objetoIndice].pontos[vt3];
    
    let p = new Ponto();
    let bari = [];
    let p_linha = new Ponto();
    for (let i = Xmin; i <= Xmax; i++) {
        p.x=i;
        p.y=Yscan;
        p.z=0;
        
        bari = acharBaricentro(p,v1,v2,v3);

        p_linha = SomaPontos(SomaPontos(multiplicarPontoPorEscalar(bari[0],vv1),multiplicarPontoPorEscalar(bari[1],vv2)),multiplicarPontoPorEscalar(bari[2],vv3));
        /*console.log(Zbuffer.length > p.x);
        console.log( Zbuffer[p.x].length > p.y);
        console.log(p_linha.z < Zbuffer[p.x][p.y]);*/
        if( Zbuffer.length > p.x && Zbuffer[p.x].length > p.y && p_linha.z < Zbuffer[p.x][p.y]){
            Zbuffer[p.y][p.x] = p_linha.z;

            let zero = new Ponto(0,0,0);
            let N = SomaVetores(SomaVetores(multiplicarVetorPorEscalar(bari[0],vv1.normal),multiplicarVetorPorEscalar(bari[1],vv2.normal)),multiplicarVetorPorEscalar(bari[2],vv3.normal));
            let V = SubtracaoPontos(zero, p_linha);
            let L = SubtracaoVetorPonto(iluminacao.Pl, p_linha);
            let R = SubtracaoVetores( multiplicarVetorPorEscalar((2*produtoInterno(L,N)),N) ,L);

            N = normalizarVetor(N);
            V = normalizarVetor(V);
            L = normalizarVetor(L);
            R = normalizarVetor(R);

            Ca = multiplicarVetorPorEscalar(iluminacao.ka, iluminacao.Ia);
            Cd = MultiplicacaoComponenteComponente(multiplicarVetorPorEscalar((iluminacao.kd * produtoInterno(N,L)),iluminacao.Od),iluminacao.Il);
            Cs = multiplicarVetorPorEscalar((iluminacao.ks*(Math.pow(produtoInterno(R,V),iluminacao.n))),iluminacao.Il);

            if(produtoInterno(N,V)<0){
                N = multiplicarVetorPorEscalar(-1,N);
            }

            if(produtoInterno(N,L)<0){
                /*iluminacao.kd = 0;
                iluminacao.ks = 0;*/
                Cd = new Vetor(0,0,0);
                Cs = new Vetor(0,0,0);
            }else{
                if(produtoInterno(R,V)<0){
                    //iluminacao.ks = 0;
                    Cs = new Vetor(0,0,0);
                }
            }
            //console.log(Ca,Cd,Cs);
            
            CIlum = SomaVetores(SomaVetores(Ca,Cd),Cs);

            if(CIlum.x > 255){
                CIlum.x = 255;
            }
            if(CIlum.y > 255){
                CIlum.y = 255;
            }
            if(CIlum.z > 255){
                CIlum.z = 255;
            }

            //console.log(CIlum);
            
            //bufferCor[p.x][p.y] = CIlum;
            bufferCor.push({x:p.x,y:p.y, cor:CIlum});
            //draw();
            //console.log(Ca,Cd,Cs);
            
        }

    }
    
}

function desenhadorBaixo(Xmin,Xmax,Yscan,vt1,vt2,v3,objetoIndice){
    let v1 = Objetos[objetoIndice].pontos_tela[vt1];
    let v2 = Objetos[objetoIndice].pontos_tela[vt2];
    //let v3 = Objetos[objetoIndice].pontos_tela[vt3];

    let vv1 = Objetos[objetoIndice].pontos[vt1];
    let vv2 = Objetos[objetoIndice].pontos[vt2];
    //let vv3 = Objetos[objetoIndice].pontos[vt3];
    let vv3 = new Ponto();
    vv3.z = v3.z;
    vv3.x = (((2*v3.x*camera.hx)/(telaX * camera.d))-(camera.hx/camera.d))*v3.z;
    vv3.y = ((camera.hx/camera.d)-((2*v3.y*camera.hy)/(telaY * camera.d)))*v3.z;
    vv3.normal = v3.normal;
    let p = new Ponto();
    let bari = [];
    let p_linha = new Ponto();
    for (let i = Xmin; i <= Xmax; i++) {
        p.x=i;
        p.y=Yscan;
        p.z=0;

        bari = acharBaricentro(p,v1,v2,v3);
        
        p_linha = SomaPontos(SomaPontos(multiplicarPontoPorEscalar(bari[0],vv1),multiplicarPontoPorEscalar(bari[1],vv2)),multiplicarPontoPorEscalar(bari[2],vv3));
        //console.log(p_linha,Zbuffer[p.x][p.y]);
        
        //console.log(Zbuffer.length > p.x);
        /*console.log(Zbuffer[p.x]);
        console.log(p_linha.z < Zbuffer[p.x][p.y]);
        console.log(p_linha.z , Zbuffer[p.x][p.y]);*/
        
        if(Zbuffer.length > p.x && Zbuffer[p.x].length > p.y && p_linha.z < Zbuffer[p.x][p.y]){
            Zbuffer[p.y][p.x] = p_linha.z;

            let zero = new Ponto(0,0,0);
            let N = SomaVetores(SomaVetores(multiplicarVetorPorEscalar(bari[0],vv1.normal),multiplicarVetorPorEscalar(bari[1],vv2.normal)),multiplicarVetorPorEscalar(bari[2],vv3.normal));
            let V = SubtracaoPontos(zero, p_linha);
            let L = SubtracaoVetorPonto(iluminacao.Pl, p_linha);
            let R = SubtracaoVetores( multiplicarVetorPorEscalar((2*produtoInterno(L,N)),N) ,L);

            N = normalizarVetor(N);
            V = normalizarVetor(V);
            L = normalizarVetor(L);
            R = normalizarVetor(R);

            Ca = multiplicarVetorPorEscalar(iluminacao.ka, iluminacao.Ia);
            Cd = MultiplicacaoComponenteComponente(multiplicarVetorPorEscalar((iluminacao.kd * produtoInterno(N,L)),iluminacao.Od),iluminacao.Il);
            Cs = multiplicarVetorPorEscalar((iluminacao.ks*(Math.pow(produtoInterno(R,V),iluminacao.n))),iluminacao.Il);

            if(produtoInterno(N,V)<0){
                N = multiplicarVetorPorEscalar(-1,N);
            }

            if(produtoInterno(N,L)<0){
                /*iluminacao.kd = 0;
                iluminacao.ks = 0;*/
                Cd = new Vetor(0,0,0);
                Cs = new Vetor(0,0,0);
            }else{
                if(produtoInterno(R,V)<0){
                    //iluminacao.ks = 0;
                    Cs = new Vetor(0,0,0);
                }
            }

            CIlum = SomaVetores(SomaVetores(Ca,Cd),Cs);
            if(CIlum.x > 255){
                CIlum.x = 255;
            }
            if(CIlum.y > 255){
                CIlum.y = 255;
            }
            if(CIlum.z > 255){
                CIlum.z = 255;
            }

            //bufferCor[p.x][p.y] = CIlum;
            bufferCor.push({x:p.x,y:p.y, cor:CIlum});
            //console.log(Ca,Cd,Cs);
            
            //draw();
        }

    }
    
}

function desenhadorCima(Xmin,Xmax,Yscan,vt1,v2,vt3,objetoIndice){
    let v1 = Objetos[objetoIndice].pontos_tela[vt1];
    //let v2 = Objetos[objetoIndice].pontos_tela[vt2];
    let v3 = Objetos[objetoIndice].pontos_tela[vt3];

    let vv1 = Objetos[objetoIndice].pontos[vt1];
    //let vv2 = Objetos[objetoIndice].pontos[vt2];
    let vv3 = Objetos[objetoIndice].pontos[vt3];
    let vv2 = new Ponto();
    vv2.z = v2.z;
    vv2.x = (((2*v2.x*camera.hx)/(telaX * camera.d))-(camera.hx/camera.d))*v2.z;
    vv2.y = ((camera.hx/camera.d)-((2*v2.y*camera.hy)/(telaY * camera.d)))*v2.z;
    vv2.normal = v2.normal;
    let p = new Ponto();
    let bari = [];
    let p_linha = new Ponto();
    for (let i = Xmin; i <= Xmax; i++) {
        p.x=i;
        p.y=Yscan;
        p.z=0;

        bari = acharBaricentro(p,v1,v2,v3);
        //console.log(vv1,vv2,vv3);
        
        p_linha = SomaPontos(SomaPontos(multiplicarPontoPorEscalar(bari[0],vv1),multiplicarPontoPorEscalar(bari[1],vv2)),multiplicarPontoPorEscalar(bari[2],vv3));
        //console.log(Zbuffer.length > p.x);
        //console.log( Zbuffer[p.x].length > p.y);
        //console.log(p_linha.z < Zbuffer[p.x][p.y]);
        
        if( Zbuffer.length > p.x && Zbuffer[p.x].length > p.y && p_linha.z < Zbuffer[p.x][p.y]){
            Zbuffer[p.y][p.x] = p_linha.z;

            let zero = new Ponto(0,0,0);
            let N = SomaVetores(SomaVetores(multiplicarVetorPorEscalar(bari[0],vv1.normal),multiplicarVetorPorEscalar(bari[1],vv2.normal)),multiplicarVetorPorEscalar(bari[2],vv3.normal));
            let V = SubtracaoPontos(zero, p_linha);
            let L = SubtracaoVetorPonto(iluminacao.Pl, p_linha);
            let R = SubtracaoVetores( multiplicarVetorPorEscalar((2*produtoInterno(L,N)),N) ,L);

            N = normalizarVetor(N);
            V = normalizarVetor(V);
            L = normalizarVetor(L);
            R = normalizarVetor(R);

            Ca = multiplicarVetorPorEscalar(iluminacao.ka, iluminacao.Ia);
            Cd = MultiplicacaoComponenteComponente(multiplicarVetorPorEscalar((iluminacao.kd * produtoInterno(N,L)),iluminacao.Od),iluminacao.Il);
            Cs = multiplicarVetorPorEscalar((iluminacao.ks*(Math.pow(produtoInterno(R,V),iluminacao.n))),iluminacao.Il);

            if(produtoInterno(N,V)<0){
                N = multiplicarVetorPorEscalar(-1,N);
            }

            if(produtoInterno(N,L)<0){
                /*iluminacao.kd = 0;
                iluminacao.ks = 0;*/
                Cd = new Vetor(0,0,0);
                Cs = new Vetor(0,0,0);
            }else{
                if(produtoInterno(R,V)<0){
                    //iluminacao.ks = 0;
                    Cs = new Vetor(0,0,0);
                }
            }

            CIlum = SomaVetores(SomaVetores(Ca,Cd),Cs);
            if(CIlum.x > 255){
                CIlum.x = 255;
            }
            if(CIlum.y > 255){
                CIlum.y = 255;
            }
            if(CIlum.z > 255){
                CIlum.z = 255;
            }

           //bufferCor[p.x][p.y] = CIlum;
           console.log(CIlum);
           bufferCor.push({x:p.x, y:p.y, cor:CIlum});
            //console.log(Ca,Cd,Cs);
            
        }

    }
    
}

function acharBaricentro(p,v1,v2,v3){
    let baricentro = [];

    let area = parseFloat(normaVetor(produtoVetorial(SubtracaoPontos(v1,v2),SubtracaoPontos(v1,v3))))/2;
    let area1 = parseFloat(normaVetor(produtoVetorial(SubtracaoPontos(p,v2),SubtracaoPontos(p,v3))))/2;
    let area2 = parseFloat(normaVetor(produtoVetorial(SubtracaoPontos(v1,p),SubtracaoPontos(v1,v3))))/2;
    let area3 = parseFloat(normaVetor(produtoVetorial(SubtracaoPontos(v1,v2),SubtracaoPontos(v1,p))))/2;
    //console.log(area, area1, area2, area3);

    let gama = (parseFloat(area3)/parseFloat(area));
    let beta = (parseFloat(area2)/parseFloat(area));
    let alfa = (parseFloat(area1)/parseFloat(area));

    baricentro[0] = alfa ;//alfa
    baricentro[1] = beta ;//beta
    baricentro[2] = gama ;//gama
   
    
    return baricentro;
}


function draw(){
    var c=document.getElementById("canvas");
    canvas.width = telaX;
    canvas.heigth = telaY;
    var ctx=c.getContext("2d");
   
    for (let i = 0; i < bufferCor.length; i++) {
        ctx.fillRect(bufferCor[i].x,bufferCor[i].y,1,1);
        ctx.fillStyle = "rgb("+bufferCor[i].cor.r+","+bufferCor[i].cor.g+","+bufferCor[i].cor.b+")";
    }
}
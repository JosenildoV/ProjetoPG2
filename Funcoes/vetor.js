class Vetor{

    constructor(x,y,z){
        this.x = x;
        this.y = y;
        this.z = z;
    }
}

function produtoVetorial(vetorA,vetorB){
    let vetorC= new Vetor();
    vetorC.x = (vetorA.y*vetorB.z) - (vetorA.z*vetorB.y);

    vetorC.y = (vetorA.z*vetorB.x) - (vetorA.x*vetorB.z);

    vetorC.z = (vetorA.x*vetorB.y) - (vetorA.y*vetorB.x);

    return vetorC;
}

function produtoInterno(vetorA,vetorB){
    let escalar;
    escalar = 0;
    escalar += (vetorA.x*vetorB.x);
    escalar += (vetorA.y*vetorB.y);
    escalar += (vetorA.z*vetorB.z);
    return escalar;
}

function MultiplicacaoComponenteComponente(vetorA,vetorB){
    let vetorC = new Vetor();

    vetorC.x = (vetorA.x*vetorB.x);
    vetorC.y = (vetorA.y*vetorB.y);
    vetorC.z = (vetorA.z*vetorB.z);
    
    return vetorC;
}

function projecaoVetores(vetorB,vetorA){
    let escalar = (produtoInterno(vetorA,vetorB)/produtoInterno(vetorB,vetorB)); 
    let proj = new Vetor();
    
    proj.x = (escalar*vetorB.x);
    proj.y = (escalar*vetorB.y);
    proj.z = (escalar*vetorB.z);

    return proj;
}

function multiplicarVetorPorEscalar(escalar, vetorA){
    let vetorB = new Vetor();
    
    vetorB.x = (escalar * vetorA.x);
    vetorB.y = (escalar * vetorA.y);
    vetorB.z = (escalar * vetorA.z);

    return vetorB;
}

function moduloVetor(vetorA){
    let modulo;
    modulo = 0;
    modulo = Math.sqrt((vetorA.x*vetorA.x)+(vetorA.y*vetorA.y)+(vetorA.z*vetorA.z));
    return modulo;
}

function normalizarVetor(vetorA){    
    let vetorB= new Vetor();
    
    let aux;
    aux = 0;
    aux = (1/moduloVetor(vetorA));
    
    vetorB = multiplicarVetorPorEscalar(aux, vetorA);

    return vetorB;
}

function OrtogonalizarVetor(vetorV, vetorN){
    let vetorOrto = new Vetor();
   
    let proj = projecaoVetores(vetorN,vetorV);
    vetorOrto.x = vetorV.x - proj.x;
    vetorOrto.y = vetorV.y - proj.y;
    vetorOrto.z = vetorV.z - proj.z;
    
    return vetorOrto;
}

function SomaVetores(vetorA,vetorB){
    let vetorC = new Vetor();
    vetorC.x = vetorA.x+vetorB.x;
    vetorC.y = vetorA.y+vetorB.y;
    vetorC.z = vetorA.z+vetorB.z;

    return vetorC;
}
function SubtracaoVetores(vetorA,vetorB){
    let vetorC = new Vetor();
    vetorC.x = vetorA.x-vetorB.x;
    vetorC.y = vetorA.y-vetorB.y;
    vetorC.z = vetorA.z-vetorB.z;

    return vetorC;
}
function SubtracaoVetorPonto(vetorA,ponto){
    let vetorC = new Vetor();
    vetorC.x = vetorA.x-ponto.x;
    vetorC.y = vetorA.y-ponto.y;
    vetorC.z = vetorA.z-ponto.z;

    return vetorC;
}

function matrizVetores(vetorA,vetorB,vetorC){
    let matriz=[];
    matriz[0]=vetorA;
    matriz[1]=vetorB;
    matriz[2]=vetorC;

    return matriz;
}

function MultMatrizVetor(matriz, vetor){
    let vetor_resultante = new Ponto();

    vetor_resultante.x = produtoInterno(matriz[0],vetor);
    vetor_resultante.y = produtoInterno(matriz[1],vetor);
    vetor_resultante.z = produtoInterno(matriz[2],vetor);

    return vetor_resultante;
}

function normaVetor(vetor){

    let norma=0;

    norma += Math.pow(vetor.x,2);
    norma += Math.pow(vetor.y,2);
    norma += Math.pow(vetor.z,2);
    norma = Math.sqrt(norma);

    return norma;
}
var arquivo_objetos = [];
var arquivo_iluminacao;
var arquivo_camera;
let objetos_carregados = false;
let camera_carregada = false;
let iluminacao_carregada = false;
var upFile = function(parametro) {

    if(parametro=="id_objetos"){
        var objectsFileSelected = document.getElementById('id_objetos').files;
        var objectsFileTobeRead = objectsFileSelected;
        let indice=0;
        var objectsFileReader = [];
        
        for (let i = 0; i < objectsFileSelected.length; i++) {
            objectsFileReader[i] = new FileReader();
        }
        
        objectsFileReader[indice].onload = function (e) {
            
            for (let j = 0; j < objectsFileSelected.length; j++) {
                
                arquivo_objetos[j] = objectsFileReader[j].result;
                let nome = objectsFileSelected[j].name;
                decomporArquivo(j,nome); 
            }
        }

        while(indice<objectsFileSelected.length){
            objectsFileReader[indice].readAsText(objectsFileTobeRead[indice]);
            indice++;
        }
    }else if(parametro=="id_camera"){
        var cameraFileSelected = document.getElementById('id_camera');
        var cameraFileTobeRead = cameraFileSelected.files[0];
        var cameraFileReader = new FileReader();

        cameraFileReader.onload = function (e){
            arquivo_camera = cameraFileReader.result;
            decomporArquivoCamera(arquivo_camera);
        }

        cameraFileReader.readAsText(cameraFileTobeRead);
    }
    else{
        var illumFileSelected = document.getElementById('id_iluminacao');
        var illumFileTobeRead = illumFileSelected.files[0];
        var illumFileReader = new FileReader();

        illumFileReader.onload = function (e){
            arquivo_iluminacao = illumFileReader.result;
            decomporArquivoIluminacao(arquivo_iluminacao);
        }

        illumFileReader.readAsText(illumFileTobeRead);
        
    }
}


function decomporArquivo(indice, nome){
    //console.log(arquivo_objetos[indice]);
    
    var linhas = arquivo_objetos[indice].split('\n');
    var qtd = linhas[0].split(' ');
    qtd[0] = parseInt(qtd[0]);
    qtd[1] = parseInt(qtd[1]);    
    
    var pts;
    var ponto;
    var pontosVet = [];
    for (let i = 1; i <= qtd[0]; i++) {
        pts = linhas[i];

        ponto = pts.split(' ');
        pontosVet[i-1] = new Ponto(parseFloat(ponto[0]),parseFloat(ponto[1]),parseFloat(ponto[2]));
    }
    
    var tri;
    var triangulo;
    var triangulosTri = [];
    for (let a = qtd[0]+1; a <= qtd[0]+qtd[1]; a++) {
        tri = linhas[a];
        
        triangulo =tri.split(' ');
        
        triangulosTri[a-1-qtd[0]] = new Triangulo(triangulo[0]-1,triangulo[1]-1,triangulo[2]-1);
        
    }

    let posicao_string = nome.indexOf(".");
    nome = nome.substring(0,posicao_string);

    Objetos[Objetos.length] = new Objeto(nome,qtd[0],qtd[1],pontosVet,triangulosTri);
    
    objetos_carregados = true;
    executar();
}

function decomporArquivoCamera(arquivo){
    let linhas = arquivo.split('\n');

    let c = linhas[0].split(' ');
    let posicao_camera = new Ponto(parseFloat(c[0]), parseFloat(c[1]), parseFloat(c[2]));

    let n = linhas[1].split(' ');    
    let vetor_N = new Vetor(parseFloat(n[0]), parseFloat(n[1]), parseFloat(n[2]));

    let v = linhas[2].split(' ');
    let vetor_V = new Vetor(parseFloat(v[0]), parseFloat(v[1]), parseFloat(v[2]));

    let outros = linhas[3].split(' ');
    let distancia = parseFloat(outros[0]);
    let hx = parseFloat(outros[1]);
    let hy = parseFloat(outros[2]);
    
    camera = new Camera(posicao_camera, vetor_N, vetor_V, distancia, hx, hy);
    camera_carregada = true;
    executar();
}

function decomporArquivoIluminacao(arquivo){
    let linhas = arquivo.split('\n');

    let l = linhas[0].split(' ');
    let posicao_luz = new Ponto(parseFloat(l[0]), parseFloat(l[1]), parseFloat(l[2]));

    let aux = linhas[1].split(' ');
    let ka = parseFloat(aux[0]);

    let corA = linhas[2].split(' ');
    let Ia = new Vetor(parseFloat(corA[0]), parseFloat(corA[1]), parseFloat(corA[2]));
    
    let aux2 = linhas[3].split(' ');
    let kd = parseFloat(aux2[0]);

    let dif = linhas[4].split(' ');
    let Od = new Vetor(parseFloat(dif[0]), parseFloat(dif[1]), parseFloat(dif[2]));

    let esp = linhas[5].split(' ');
    let ks = parseFloat(esp[0]);
    
    let corL = linhas[6].split(' ');
    let Il = new Vetor(parseFloat(corL[0]), parseFloat(corL[1]), parseFloat(corL[2]));

    let rugo = linhas[7].split(' ');
    let n = parseFloat(rugo[0]);

    iluminacao = new Iluminacao(posicao_luz, ka, Ia, kd, Od, ks, Il, n);
    
    iluminacao_carregada = true;
    executar();
}

function executar(){
    if(objetos_carregados && camera_carregada && iluminacao_carregada){
        Inicializar();
    }
}
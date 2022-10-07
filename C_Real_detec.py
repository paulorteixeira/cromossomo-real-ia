import numpy as np
import random as rd
import matplotlib.pyplot as plt

##definindo os limites das variaveis
lim_sup_var1  = 12
lim_inf_var1  = -10

lim_sup_var2  = 12
lim_inf_var2  = -10


rd.seed(7)


#definindo a função a ser minimizada
#def f (x,y):
#    return (15+((x-3)**2)/2+((y-3)**2)/2-2*(np.sin(4*x-3)+np.sin(4*y-3)))

#definindo a função a ser minimizada
def f (x,y):
    return (150-x**2-y**2)

    
 #   real(0.95,0.5,     100,      200,          True,     2,  15,    True)
def real(pc,   pm,geracoes,população,roleta_torneio,pontos,kpop, elitismo):
    tamCromossomo = 2 ##tamanho do cromossomo -> n vars


    pc = pc          ##prop de cruzamento
    pm = pm          ##prop de mutações
    numgeracoes  = geracoes ## numero de gerações
    tamPopulacao = população  ## tamanho da populacao
    roleta_torneio = roleta_torneio
    pontosCruzamento = pontos
    kpop = kpop
    elitismo = elitismo


    ## geração aleatoria da populaçao inicial

    p = np.zeros((tamPopulacao,tamCromossomo))##populaçao
    for i in range(tamPopulacao):
        for j in range (tamCromossomo):
            a = rd.uniform(0,1)
            
            if(j==0):
                p[i][j] = lim_inf_var1 + a * (lim_sup_var1-lim_inf_var1)
            if(j==1):
                p[i][j] = lim_inf_var2 + a * (lim_sup_var2-lim_inf_var2)
            

    ## Criação de vars do AG
    ind = np.zeros(tamCromossomo)  
    individuo_var1 = np.zeros(tamPopulacao) ##valores reais
    individuo_var2 = np.zeros(tamPopulacao) ##valores reais

    Aptidao = np.zeros(tamPopulacao)
    novageracao = np.zeros((tamPopulacao,tamCromossomo))
    geracoes = 0

    ## iniciando o AG
    while (geracoes<=numgeracoes):
        novosindividuos = 0
        while (novosindividuos<(tamPopulacao-1)):
            ##Tranformando individuos de bin para real
            for i in range(tamPopulacao):
                ind[:] = p[i,:]
                conv = 0
                for j in range(tamCromossomo):
                    if(j==0):
                        individuo_var1[i] = ind[j]
                    if(j==1):
                        individuo_var2[i] = ind[j]
                        
            ##calculo da aptidao dos individuos
            TotalAptidao = 0
            for i in range(tamPopulacao):
                Aptidao[i] = 1/f( individuo_var1[i], individuo_var2[i])## adaptação da função f(x)
                TotalAptidao = Aptidao[i]+TotalAptidao
                
            ## Selecao dos pais para cruzamento  
            ## identificar a probabilidade de cada individuo 
            pic = np.zeros(tamPopulacao)   
            pitotal = np.zeros(tamPopulacao)
            pic = (1/TotalAptidao)*Aptidao
            #criando roleta
            if(roleta_torneio == True):
                for i in range(tamPopulacao):
                    if(i == 0):
                        pitotal[i] = pic[i]
                    else:
                        pitotal[i] = pic[i]+pitotal[i-1]
                ## sorteando os pais de acordo com a prob
                roleta1 = rd.uniform(0,1)
                i = 0
                while (roleta1>pitotal[i]):
                    i = i+1
                pai1 = i

                roleta2 = rd.uniform(0,1)
                i = 0
                while (roleta2>pitotal[i]):
                    i = i+1
                pai2 = i

                while(pai2 == pai1):
                    roleta2 = rd.uniform(0,1)
                    i = 0
                    while (roleta2>pitotal[i]):
                        i = i+1
                    pai2 = i
            ##fim roleta
            else:
                auxii = 0
                maxApt1 = 0
                maxApt2 = 0
                k_sorteado = np.zeros(kpop)
                i = 0
                while(i< kpop): 
                    k_sorteado[i] = int(round(1+(tamCromossomo-2)*rd.uniform(0,1)))
                    i = i +1
               
                i=0
                while(i< kpop):
                    
                    opi = int(k_sorteado[i])
                    if(maxApt1<Aptidao[opi]):
                        maxApt1 = Aptidao[opi]
                        pai1 = int(k_sorteado[i])
                    i = i +1
                    
                i = 0
                while(i< kpop):
                    opi = int(k_sorteado[i])
                    if(maxApt2<Aptidao[opi] and maxApt2!=maxApt1):
                        maxApt2 = Aptidao[opi]
                        pai2 = int(k_sorteado[i])
                    i = i +1  

                
            ## RADCLIFF
            if(pontosCruzamento == 1):
                if(pc>rd.uniform(0,1)):
                
                    geneF1 = np.zeros(tamCromossomo)
                    geneF2 = np.zeros(tamCromossomo)
                    
                    for i in range(tamCromossomo):
                        beta = rd.uniform(0,1)
                        geneF1[i] = beta*p[pai1][i] + (1-beta)*p[pai2][i]
                        geneF2[i] = (1-beta)*p[pai1][i] + beta*p[pai2][i]              
        
                    novageracao[novosindividuos,:] = geneF1
                    novosindividuos = novosindividuos+1    
                    novageracao[novosindividuos,:] = geneF2                        
                    novosindividuos = novosindividuos+1
            else:
                ##wright
                if(pc>rd.uniform(0,1)):
                    a=b=c=-1024
                    geneFilho1 = np.zeros(tamCromossomo)
                    geneFilho2 = np.zeros(tamCromossomo)
                    geneF1 = np.zeros(tamCromossomo)
                    geneF2 = np.zeros(tamCromossomo)
                    geneF3 = np.zeros(tamCromossomo)
                    
                    for i in range(tamCromossomo):
                        if(i==0):
                            
                            geneF1[i] =  0.5*p[pai1][i] + 0.5*p[pai2][i]
                                
                            geneF2[i] =  1.5*p[pai1][i] - 0.5*p[pai2][i]
                                                        
                            geneF3[i] = -0.5*p[pai1][i] + 1.5*p[pai2][i]
                                                        
                        if(i==1):
                            
                            geneF1[i] =  0.5*p[pai1][i] + 0.5*p[pai2][i]
                                
                            geneF2[i] =  1.5*p[pai1][i] - 0.5*p[pai2][i]
                                
                            geneF3[i] = -0.5*p[pai1][i] + 1.5*p[pai2][i]
                            
                    valido1 =valido2=valido3 = False
                    if( geneF1[0] >= lim_inf_var1 and geneF1[0] <= lim_sup_var1 and geneF1[1] >= lim_inf_var2 and geneF1[1] <= lim_sup_var2):
                         valido1 = True
                         a = 1/f(geneF1[0],geneF1[1])
                    if( geneF2[0] >= lim_inf_var1 and geneF2[0] <= lim_sup_var1 and geneF2[1] >= lim_inf_var2 and geneF2[1] <= lim_sup_var2):
                         valido2 = True
                         b = 1/f(geneF2[0],geneF2[1])
                    if( geneF3[0] >= lim_inf_var1 and geneF3[0] <= lim_sup_var1 and geneF3[1] >= lim_inf_var2 and geneF3[1] <= lim_sup_var2):
                         valido3 = True
                         c = 1/f(geneF3[0],geneF3[1])
                         
                    arr = [a,b,c]
                        
                    arr1 = np.sort(arr)
                    
                    if((not valido1 and not valido2 and not valido3)):
                        for i in range(tamCromossomo):
                            beta = rd.uniform(0,1)
                            geneFilho1[i] = beta*p[pai1][i] + (1-beta)*p[pai2][i]
                            geneFilho2[i] = (1-beta)*p[pai1][i] + beta*p[pai2][i]
                        
                    else:
                        if(arr1[2] == a and valido1):
                             geneFilho1 = geneF1
                        if(arr1[2] == b and valido2):
                             geneFilho1 = geneF2
                        if(arr1[2] == c and valido3):
                             geneFilho1 = geneF3
                        
                        if(arr1[1] == a and valido1):
                            geneFilho2 = geneF1
                        if(arr1[1] == b and valido2):
                            geneFilho2 = geneF2
                        if(arr1[1] == c and valido3):
                            geneFilho2 = geneF3
                
                    novageracao[novosindividuos,:] = geneFilho1                  
                    novosindividuos = novosindividuos+1
                    novageracao[novosindividuos,:] = geneFilho2
                    novosindividuos = novosindividuos+1

                
            ##mutação
            if(pm>rd.uniform(0,1)):
                d = round(1+(tamCromossomo-2)*rd.uniform(0,1)) 
                for j in range (tamCromossomo):
                    a = rd.uniform(0,1)
                    
                    if(j==0):
                        novageracao[novosindividuos-2][d] = lim_inf_var1 + a * (lim_sup_var1-lim_inf_var1)
                    if(j==1):
                        novageracao[novosindividuos-1][d] = lim_inf_var2 + a * (lim_sup_var2-lim_inf_var2)
                


        indice = Aptidao.argmax()
        elem_var1 = individuo_var1[indice]
        elem_var2 = individuo_var2[indice]

        if(elitismo == True):
            novageracao[0] = p[indice]
        
        if(geracoes>0):
            if(elem_var1>melhoraptidao_var1):
                melhoraptidao_var1 = elem_var1
                indicemelhoraptidao = indice
                melhor_var1 = individuo_var1[indice]
                #y = -abs(melhorx*np.sin(np.sqrt(abs(melhorx))))
                melhorgeracao = geracoes
                populacao = novageracao

                
            if(elem_var2>melhoraptidao_var2):
                melhoraptidao_var2 = elem_var2
                indicemelhoraptidao = indice
                melhor_var2 = individuo_var2[indice]
                #y = -abs(melhorx*np.sin(np.sqrt(abs(melhorx))))
                melhorgeracao = geracoes
                populacao = novageracao
        else:
            melhoraptidao_var1 = elem_var1
            melhoraptidao_var2 = elem_var2

            indicemelhoraptidao = indice
            melhor_var1 = individuo_var1[indice]
            melhor_var2 = individuo_var2[indice]

            zeta = f(melhor_var1,melhor_var2)
            #y = -abs(melhorx*np.sin(np.sqrt(abs(melhorx))))
            melhorgeracao = 0
            populacao = p
        
            
        p = novageracao                        
        geracoes = geracoes + 1
        print('melhor geracao:',melhorgeracao)
        print('melhor x:', melhor_var1)
        print('melhor y:',melhor_var2)

        print('f(x,y):',f(melhor_var1,melhor_var2)  )


        
    #gerando o grafico //inf- sup-passos
    vetorx = np.linspace(0,450,100)
    vetory = -abs(vetorx*np.sin(np.sqrt(abs(vetorx))))
    plt.plot(vetorx,vetory)
    plt.plot(melhor_var1,melhor_var2,'o',color='red')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()


    print('melhor geracao:', melhorgeracao)
    print('melhor x:', melhor_var1)
    print('melhor y:', melhor_var2)

    print('f(x,y):', f(melhor_var1,melhor_var2)  )

    return 'melhor geracao:'+ str(melhorgeracao)+ '\nmelhor x:'+ str(melhor_var1)+'\nmelhor y:'+ str(melhor_var2)+ '\nf(x,y):'+str(f(melhor_var1,melhor_var2))      
##real(pc,pm,geracoes,população,roleta_torneio,pontos,kpop, elitismo)
print('calculando')
real(0.95,0.5,2,10,True,2,15, True)

from django.contrib import admin
from django.urls import path, include


admin.site.site_header = "Панель управления"

urlpatterns = [
    path('admin/', admin.site.urls),
    path('main/', include('main.urls')),
    path('accounts/',include('django.contrib.auth.urls')),
    ]

